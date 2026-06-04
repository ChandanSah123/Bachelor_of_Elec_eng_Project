
import os
import re
import numpy as np
import scipy.io as sio

# Attempt to import PSS/E modules (edit if your module names differ)
import psse3603 as psse_mod    # adjust to your installed psse module name if needed
import psspy
import dyntools

# ----------------- USER CONFIG -----------------
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC programs"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")   # your .raw
dyr_file = os.path.join(work_dir, "ieee9bus.dyr")   # optional .dyr (leave "" if none)
txt_file = os.path.join(work_dir, "Ybus_Export.txt")
mat_file = os.path.join(work_dir, "Y_all.mat")

# Fault & clearing specification (bus numbers use PSS/E bus numbering from the raw file)
fault_bus = 7          # bus where 3-phase fault is applied (for our numeric Ybuild we will remove this bus to represent solid ground)
trip_line_from = 7     # line 'from' bus to trip when clearing (for the post-fault)
trip_line_to   = 5     # line 'to' bus to trip when clearing

# Optional Kron reduction parameters:
# If you want Kron reduction to internal generator nodes behind Xd', specify:
# gen_buses: list of generator bus numbers (in same numbering as .raw)
# xd_prime: list or array with corresponding transient reactances Xd' (in pu, same system base as case)
# If left empty, Kron reduction step is skipped.
gen_buses = []   # e.g. [1,2,3]
xd_prime  = []   # e.g. [0.146, 0.8958, 1.3125]

# ------------------------------------------------

def init_psse_and_read_case():
    # init PSS/E
    psspy.psseinit(200)  # 200 buses workspace (adjust if needed)
    # read raw
    ierr = psspy.read(0, raw_file)
    if ierr != 0:
        raise RuntimeError(f"psspy.read returned error {ierr} for {raw_file}")
    # optional dynamic file
    if dyr_file:
        psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

def solve_and_convert():
    # Solve power flow (to get voltages for conversion)
    ierr_fnsl = psspy.fnsl([0,0,0,1,1,0,99,0])
    if isinstance(ierr_fnsl, tuple):
        ierr_fnsl = ierr_fnsl[0]
    if ierr_fnsl > 1:
        print("Warning: fnsl returned", ierr_fnsl)
    # Convert generators and loads to Norton/admittance representations (CONG/CONL)
    psspy.cong(0)
    # Convert loads in 3 steps to get constant admittance representation (100% for P and Q -> admittance)
    psspy.conl(0,1,1,[0,0],[100.0,0.0,0.0,100.0])
    psspy.conl(0,1,2,[0,0],[100.0,0.0,0.0,100.0])
    psspy.conl(0,1,3,[0,0],[100.0,0.0,0.0,100.0])
    # Order / factorize / type (needed before output_y_matrix sometimes)
    # Note: some PSS/E installations may block or require a specific call-order — errors do not always mean failure.
    try:
        psspy.ordr(0)
        psspy.fact()
        psspy.tysl(0)
    except Exception:
        # If any of these fail due to CONL multi-step messages, it's usually safe to continue if output_y_matrix works.
        pass

def export_ybus_to_text():
    # Use PSS/E API to dump Y-matrix to a text file.
    # Notes: output_y_matrix arguments: (sid, all, ties, out, ofile)
    # out=0 -> write to file; ofile required.
    ierr = psspy.output_y_matrix(0, 1, 0, 0, txt_file)
    if isinstance(ierr, tuple):
        ierr = ierr[0]
    if ierr != 0:
        raise RuntimeError(f"output_y_matrix returned error code {ierr}. Check PSS/E permissions and args.")
    if not os.path.exists(txt_file):
        raise FileNotFoundError(f"Expected exported text file not found: {txt_file}")

def parse_ybus_text_to_numpy():
    # Parse the text file produced by output_y_matrix, expecting lines with:
    # i j real imag
    Yreal = {}
    Yimag = {}
    max_bus = 0
    float_re = re.compile(r"[-+]?\d*\.\d+|\d+")
    with open(txt_file, "r") as f:
        for line in f:
            parts = float_re.findall(line)
            if len(parts) >= 4:
                try:
                    i = int(parts[0])
                    j = int(parts[1])
                    real = float(parts[2])
                    imag = float(parts[3])
                except Exception:
                    continue
                max_bus = max(max_bus, i, j)
                Yreal[(i, j)] = real
                Yimag[(i, j)] = imag
    # build symmetric full matrix (1-based -> 0-based)
    Y = np.zeros((max_bus, max_bus), dtype=np.complex128)
    for (i, j), rv in Yreal.items():
        iv = Yimag[(i, j)]
        Y[i-1, j-1] = rv + 1j*iv
        # mirror (in case file only contains upper/lower)
        Y[j-1, i-1] = rv + 1j*iv
    return Y

def build_faulted_physical(Y_full, fault_bus_num):
    idx = fault_bus_num - 1
    Yf = np.delete(Y_full, idx, axis=0)
    Yf = np.delete(Yf, idx, axis=1)
    return Yf

def build_postfault_physical_remove_line(Y_full, a_bus, b_bus):
    idx_a = a_bus - 1
    idx_b = b_bus - 1
    Yp = Y_full.copy()
    # If there is already zero entry, warn
    y_off = Yp[idx_a, idx_b]
    if abs(y_off) < 1e-12:
        # already 0 — nothing to remove
        return Yp
    # compute line admittance (note: offdiag = -y_line)
    y_line = -Yp[idx_a, idx_b]
    # subtract y_line from both diagonals
    Yp[idx_a, idx_a] = Yp[idx_a, idx_a] - y_line
    Yp[idx_b, idx_b] = Yp[idx_b, idx_b] - y_line
    # zero the off-diagonals between a and b
    Yp[idx_a, idx_b] = 0.0 + 0.0j
    Yp[idx_b, idx_a] = 0.0 + 0.0j
    return Yp

def kron_reduce_internal(Y_physical, gen_buses_list, xd_prime_list):
    K = len(gen_buses_list)
    N = Y_physical.shape[0]
    if K == 0:
        raise ValueError("No generator buses provided for Kron reduction.")
    if K != len(xd_prime_list):
        raise ValueError("Generator buses and xd_prime lengths mismatch.")
    # For generality we expect that user supplies positions of generator terminal buses (1-based indices relative to Y_physical).
    # But often Y_physical corresponds to system buses 1..N. Map those bus numbers -> indices (0-based).
    gen_idx = [b-1 for b in gen_buses_list]
    # Build Y_pre-like augmented matrix blocks (using complex admittance 1/(j*xd') = -1j/xd')
    # YA (KxK) = diag(1/(j*xd'))
    YA = np.zeros((K, K), dtype=np.complex128)
    for i, xd in enumerate(xd_prime_list):
        YA[i, i] = 1.0/(1j*xd)  # = -1j/xd
    # YB (K x N): coupling between internal nodes and terminal buses
    YB = np.zeros((K, N), dtype=np.complex128)
    for i, idx in enumerate(gen_idx):
        YB[i, idx] = -YA[i, i]   # -1/(j*xd)
    # YC = YB^T (N x K)
    YC = YB.T
    # YD = Y_physical (N x N)
    YD = Y_physical
    # Solve for YD\YC
    try:
        YD_inv_YC = np.linalg.solve(YD, YC)   # N x K
    except np.linalg.LinAlgError as e:
        raise RuntimeError("Kron reduction failed: Y_physical singular") from e
    Yint = YA - (YB @ YD_inv_YC)
    return Yint

# === Main execution ===
if __name__ == "__main__":
    # sanity checks
    if not os.path.exists(raw_file):
        raise FileNotFoundError(f"raw file not found: {raw_file}")
    # init and read
    print("Initializing PSS/E and reading case...")
    init_psse_and_read_case()

    # solve and convert (to get admittances)
    print("Solving power flow and converting network...")
    solve_and_convert()

    # export ybus to text and parse
    print("Exporting PRE-FAULT Ybus to text...")
    export_ybus_to_text()
    print("Parsing text -> numeric Ybus...")
    Y_pre = parse_ybus_text_to_numpy()
    print(f"Y_pre shape: {Y_pre.shape}")

    # build faulted physical by removing the fault bus
    print(f"Building FAULTED Y by removing bus {fault_bus} ...")
    Y_fault = build_faulted_physical(Y_pre, fault_bus)
    print(f"Y_fault shape: {Y_fault.shape}")

    # build post-fault physical by tripping specified line (a,b)
    print(f"Building POST-FAULT Y by tripping line {trip_line_from}-{trip_line_to} ...")
    Y_post = build_postfault_physical_remove_line(Y_pre, trip_line_from, trip_line_to)
    print(f"Y_post shape: {Y_post.shape}")

    # Optional Kron reductions (if requested)
    do_kron = (len(gen_buses) > 0 and len(xd_prime) > 0)
    if do_kron:
        print("Performing Kron reduction to internal generator nodes (pre, faulted, post)...")
        # IMPORTANT: the gen_buses entries here must refer to the *bus indices in Y_pre*
        try:
            Yint_pre = kron_reduce_internal(Y_pre, gen_buses, xd_prime)
        except Exception as e:
            print("Kron reduction pre-fault failed:", e)
            Yint_pre = None
        # For faulted Y we removed the fault bus -> if gen_buses include that bus you must adjust indices.
        # For simplicity we will attempt to map generator bus numbers to indices in the reduced matrix (skip if mapping invalid)
        try:
            # Build a physical Y for fault by removing the bus; then Kron reduce using same gen list except the faulted bus if present
            gen_buses_fault = [b for b in gen_buses if b != fault_bus]
            xd_fault = [xd_prime[i] for i,b in enumerate(gen_buses) if b != fault_bus]
            if len(gen_buses_fault) > 0:
                Yint_fault = kron_reduce_internal(Y_fault, gen_buses_fault, xd_fault)
            else:
                Yint_fault = None
        except Exception as e:
            print("Kron reduction for faulted network failed:", e)
            Yint_fault = None

        # Post-fault (line removed). If gen buses unchanged we can reduce Y_post directly.
        try:
            Yint_post = kron_reduce_internal(Y_post, gen_buses, xd_prime)
        except Exception as e:
            print("Kron reduction post-fault failed:", e)
            Yint_post = None
    else:
        Yint_pre = Yint_fault = Yint_post = None

    # Save results to mat
    print(f"Saving all matrices to MATLAB file: {mat_file} ...")
    mdict = {
        "Y_pre": Y_pre,
        "Y_fault": Y_fault,
        "Y_post": Y_post
    }
    if Yint_pre is not None:
        mdict["Yint_pre"] = Yint_pre
    if Yint_fault is not None:
        mdict["Yint_fault"] = Yint_fault
    if Yint_post is not None:
        mdict["Yint_post"] = Yint_post

    sio.savemat(mat_file, mdict)
    print("SUCCESS: All Y-bus matrices saved.")
