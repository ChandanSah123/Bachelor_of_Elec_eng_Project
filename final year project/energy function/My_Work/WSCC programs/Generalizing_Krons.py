# build_reduced_Y.py
# Generalized script: read .raw, build Ybus, detect generators, add internal nodes behind Xd',
# compute Kron-reduced Y_int for pre-fault, faulted (bus grounded), and post-fault (line trip).
# Saves results to a .mat file.

import os
import re
import sys
import numpy as np
import scipy.io as sio

# try to import psspy
try:
    import psse3603 as _psse_mod  # optional local name
except Exception:
    pass
try:
    import psspy
except Exception as e:
    raise RuntimeError("psspy import failed. Run this script in PSS/E python environment.") from e

# ------------------- user config -------------------
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC programs"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")   # change if needed
dyr_file = os.path.join(work_dir, "ieee9bus.dyr")   # optional (not required)
txt_ybus = os.path.join(work_dir, "Ybus_temp.txt")
out_mat = os.path.join(work_dir, "Reduced_Y_all.mat")

# Fault & post-fault configuration (example)
faulted_bus_user = None   # set to integer bus number to ground a bus; None => no additional fault
line_trip = None          # set to tuple (from_bus, to_bus) to remove that line after clearing, e.g. (5,7)

# If automatic Xd' extraction fails, supply transient reactance vector here (index order will be generator order detected below)
Xd_primes_manual = None   # e.g. [0.25, 0.3, 0.35]  OR None to rely on automatic detection or fallback default

# default fallback transient reactance (pu) if not found
DEFAULT_XD_PRIME = 0.3

# ------------------- helper functions -------------------
def safe_call_psse(fn, *args, **kwargs):
    """Call a psspy function and return (ierr, result) where possible; for functions returning int code just return code."""
    try:
        out = fn(*args, **kwargs)
        return out
    except Exception as e:
        # many psspy functions either return (ierr, val) or raise; normalize
        raise

def parse_output_y_file(txt_path):
    """Parse the PSS/E output_y_matrix text file into a full dense Ybus (complex numpy).
       Returns Y (N x N), and mapping of bus numbers present (list).
    """
    if not os.path.isfile(txt_path):
        raise FileNotFoundError(f"Ybus text file not found: {txt_path}")

    entries = []  # tuples (i, j, real, imag)
    with open(txt_path, 'r') as f:
        for line in f:
            # each data line typically: "    i,    j, real, imag"
            # capture numbers (integers and floats, positive/negative)
            parts = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            if len(parts) >= 4:
                try:
                    i = int(parts[0])
                    j = int(parts[1])
                    real = float(parts[2])
                    imag = float(parts[3])
                    entries.append((i, j, real, imag))
                except:
                    continue

    if not entries:
        raise RuntimeError("No numerical entries parsed from Ybus text file.")

    maxbus = max(max(e[0], e[1]) for e in entries)
    Y = np.zeros((maxbus, maxbus), dtype=complex)
    bus_present = set()
    for i, j, r, im in entries:
        Y[i-1, j-1] = r + 1j*im
        Y[j-1, i-1] = r + 1j*im
        bus_present.add(i); bus_present.add(j)
    return Y, sorted(bus_present)

def try_extract_generators_from_psse():
    """Try to use psspy API to list generators and extract Xd' (if available).
       Returns list of generator bus numbers and list of Xd' values (or None).
    """
    gen_buses = []
    xdp_list = []
    # Common approach: use psspy.macdat or psspy.machine_data_2 or psspy.machdat
    # We'll try a few calls in try/except blocks; if they fail we fall back to raw parsing.
    try:
        # attempt to get list of machines using psspy._i and psspy.getmachine... (psspy API differs)
        # psspy.machine_data_2 returns a tuple (ierr, [fieldlist]) in some versions. We attempt to loop over bus numbers.
        ierr, buses = psspy.abusint(-1, 'NUMBER')  # try to get list of buses: this may raise depending on version
    except Exception:
        buses = None

    # Try a safe search: query generator count using psspy.amachcount or psspy.agenbus / psspy.mach... Not portable;
    # Instead try iterating buses and asking if machine exists at that bus using psspy.macdat
    try:
        # some PSS/E versions expose psspy.macdat(busnum, 'parameter') but this is version dependent.
        # We'll attempt a robust loop over bus numbers 1..maxbus and check for a machine by requesting generator real power
        ierr, bus_list = psspy.abusint(-1, 'NUMBER')  # maybe returns list of bus numbers
        if isinstance(bus_list, (list, tuple, np.ndarray)):
            maxbus = max(bus_list)
        else:
            maxbus = 0
    except Exception:
        # fallback: small search
        maxbus = 500

    for b in range(1, maxbus+1):
        try:
            # try to get list of machines at bus b: use psspy.machine_busnumber or psspy.machines? API differs
            # We'll attempt psspy.machdat which exists in many versions
            # signature (ierr, data) = psspy.machdat(busnum, machine_id, 'parameter')
            # We'll check if there is machine id '1' at bus b
            ierr, _ = psspy.machdat(b, '1', ['status'])  # 'status' commonly available
            if ierr == 0:
                # machine found
                gen_buses.append(b)
                # try to get Xd' (common parameter names: 'xd_p', 'xd' or index based)
                # attempt 'xd_p' parameter name
                try:
                    ierr2, xdpv = psspy.machdat(b, '1', ['xd_p'])
                    if ierr2 == 0 and xdpv is not None:
                        # machdat often returns a list of lists so normalize
                        val = xdpv[0] if isinstance(xdpv, (list,tuple)) else xdpv
                        xdp_list.append(float(val))
                        continue
                except Exception:
                    pass
                # try other parameter labels or API calls if needed (skip for brevity)
                xdp_list.append(None)
        except Exception:
            # not found or API not available
            pass

    if not gen_buses:
        return None, None
    # normalize xdp_list: if all None return None
    if all(v is None for v in xdp_list):
        return gen_buses, None
    return gen_buses, xdp_list

def parse_raw_for_generators(raw_path):
    """Fallback parser: read RAW file, find MACHINE or GENERATOR data blocks and extract bus numbers and attempt Xd'.
       RAW file formats vary; this function tries common patterns:
       - Searches for 'MACHINE DATA FOLLOWS' or 'GENERATOR DATA FOLLOWS'
       - Then reads lines under that heading until an empty line or 'END'
       - Each machine line typically begins with bus number.
    """
    gens = []
    xdp = []
    with open(raw_path, 'r') as f:
        txt = f.read().splitlines()
    start_idx = None
    for i, line in enumerate(txt):
        if re.search(r"MACHINE\s+DATA", line, re.IGNORECASE) or re.search(r"GENERATOR\s+DATA", line, re.IGNORECASE):
            start_idx = i+1
            break
    if start_idx is None:
        # try another pattern: RAW v33 has a 'MACHINE DATA FOLLOWS' then one-line headers
        for i, line in enumerate(txt):
            if line.strip().upper().startswith('MACHINE DATA'):
                start_idx = i+1
                break
    if start_idx is None:
        # no marker found
        return [], []

    # read subsequent lines until blank or '0' alone, or next section (lines starting with '0 /')
    for j in range(start_idx, len(txt)):
        line = txt[j].strip()
        if not line:
            break
        # some raw formats use commas and continued lines, we attempt to parse first token as bus number
        # skip header-like lines
        if line.startswith('0') and ('/' in line or line.strip() == '0'):
            break
        tokens = re.split(r'[,\s]+', line)
        if len(tokens) < 1:
            continue
        try:
            busnum = int(tokens[0])
            gens.append(busnum)
            # attempt to find an xdp in the tokens (scan tokens for a number in plausible Xd' range)
            found_xd = None
            for tok in tokens[1:]:
                try:
                    val = float(tok)
                    if 0.0 < val < 10.0:  # plausible pu reactance/h parameter — heuristic
                        found_xd = val
                        break
                except:
                    continue
            xdp.append(found_xd)
        except:
            # not a machine data line
            continue
    return gens, xdp

def kron_reduction(Y_aug, keep_indices):
    """Perform Kron reduction on Y_aug, retaining rows/cols in keep_indices (list of indices 0-based),
       i.e. eliminate all nodes not in keep_indices. Returns reduced Y matrix sized len(keep_indices).
       Equivalent to Y_red = Y_kk - Y_krest * (Y_restrest^{-1}) * Y_restk
    """
    N = Y_aug.shape[0]
    keep = np.array(keep_indices, dtype=int)
    rest = np.array([i for i in range(N) if i not in keep], dtype=int)
    if rest.size == 0:
        return Y_aug[np.ix_(keep, keep)].copy()
    Y_kk = Y_aug[np.ix_(keep, keep)]
    Y_krest = Y_aug[np.ix_(keep, rest)]
    Y_restk = Y_aug[np.ix_(rest, keep)]
    Y_restrest = Y_aug[np.ix_(rest, rest)]
    # solve rather than invert
    try:
        sol = np.linalg.solve(Y_restrest, Y_restk)
        Yred = Y_kk - Y_krest.dot(sol)
        return Yred
    except np.linalg.LinAlgError:
        raise RuntimeError("Kron reduction failed: rest-rest matrix singular.")

# ------------------- main sequence -------------------
if __name__ == "__main__":
    if not os.path.isfile(raw_file):
        raise FileNotFoundError(f"Raw file not found: {raw_file}")

    # init psspy
    psspy.psseinit(50)

    # read raw (0 = default case id)
    ierr = psspy.read(0, raw_file)
    if ierr != 0:
        print("Warning: psspy.read returned code", ierr)

    # solve power flow
    ierr_fnsl = psspy.fnsl([0,0,0,1,1,0,99,0])
    if isinstance(ierr_fnsl, tuple):  # some versions return (ierr, ...)
        ierr_code = ierr_fnsl[0]
    else:
        ierr_code = ierr_fnsl
    print("fnsl returned", ierr_code)
    if ierr_code != 0:
        print("Warning: fnsl returned non-zero. Results may be invalid.")

    # convert network for admittance extraction
    psspy.cong(0)
    # convert loads to admittances (3-step recommended)
    psspy.conl(0,1,1,[0,0],[100.0,0.0,0.0,100.0])
    psspy.conl(0,1,2,[0,0],[100.0,0.0,0.0,100.0])
    psspy.conl(0,1,3,[0,0],[100.0,0.0,0.0,100.0])

    # attempt to output Y-bus to text file
    # please note: output_y_matrix signature differs by version; using (sid, all, ties, out, ofile)
    print("Writing Y-bus text file:", txt_ybus)
    ierr_out = psspy.output_y_matrix(0, 1, 0, 0, txt_ybus)  # out=0 -> write to file OFILE
    if isinstance(ierr_out, tuple):
        ierr_out = ierr_out[0]
    print("output_y_matrix returned", ierr_out)
    if ierr_out != 0:
        raise RuntimeError(f"output_y_matrix failed with code {ierr_out}. Check API args for this PSS/E version.")

    # parse the Y-bus text file
    Ybus_full, bus_list = parse_output_y_file(txt_ybus)
    print(f"Parsed Y-bus full matrix size: {Ybus_full.shape}, buses found: {len(bus_list)}")

    # Try to get generator buses and Xd' using PSS/E API
    gen_buses, xdp_from_psse = try_extract_generators_from_psse()
    if not gen_buses:
        print("Could not extract generators via PSS/E API; falling back to parsing raw file.")
        gen_buses, xdp_from_raw = parse_raw_for_generators(raw_file)
        if gen_buses:
            print("Generators detected from raw file:", gen_buses)
        else:
            raise RuntimeError("No generators detected automatically. Please provide generator bus list manually.")
        # merge xdp info
        if any(v is not None for v in xdp_from_raw):
            xdp_list = xdp_from_raw
        else:
            xdp_list = None
    else:
        print("Generators detected via PSS/E API:", gen_buses)
        xdp_list = xdp_from_psse

    # If user provided manual Xd' vector, use it
    if Xd_primes_manual is not None:
        if len(Xd_primes_manual) != len(gen_buses):
            print("Warning: manual Xd' length mismatch with number of detected generators. Ignoring manual Xd'.")
        else:
            xdp_list = list(Xd_primes_manual)

    # Build Xd' vector (fallback to default where missing)
    if xdp_list is None:
        xdp_list = [None] * len(gen_buses)

    for i in range(len(xdp_list)):
        if xdp_list[i] is None:
            print(f"Warning: Xd' for generator at bus {gen_buses[i]} not found; using default {DEFAULT_XD_PRIME}.")
            xdp_list[i] = DEFAULT_XD_PRIME

    # Now construct augmented admittance matrix that includes internal nodes for each machine
    # Approach:
    # - Ybus_full is NxN for the full network (bus numbering 1..N)
    # - We'll create internal nodes: for M generators we will add M internal nodes (one per gen)
    # - For each generator at terminal bus b (1-based), we:
    #   - subtract nothing from Ybus_full (terminal bus stays)
    #   - add connection between internal node and terminal bus with admittance 1/(j*Xd')
    #
    # Final augmented ordering: [internal_nodes (M)] then [network buses 1..N]  -> size (M+N)x(M+N)
    M = len(gen_buses)
    N = Ybus_full.shape[0]
    print(f"Detected {M} machines; building augmented Y of size {M+N}")

    Y_aug = np.zeros((M+N, M+N), dtype=complex)

    # put network Y into bottom-right
    Y_aug[M:M+N, M:M+N] = Ybus_full.copy()

    # add internal node self-admittances and coupling
    for idx, busnum in enumerate(gen_buses):
        terminal_idx = busnum - 1  # zero-based index in Ybus_full
        xd_p = float(xdp_list[idx])
        # admittance between internal node and terminal bus (j*Xd' in denominator)
        # transient reactance Xd' is purely imaginary j*Xd', so admittance = 1/(j*Xd') = -1j / Xd'
        y_adm = -1j / xd_p
        # internal node index = idx (0..M-1)
        int_i = idx
        term_j = M + terminal_idx
        # fill submatrices:
        Y_aug[int_i, int_i] += y_adm      # internal self admittance
        Y_aug[term_j, term_j] += y_adm    # add coupling to terminal diagonal as well (since admittance appears in the network)
        Y_aug[int_i, term_j] += -y_adm    # off-diagonal (internal -> terminal) should be -y_adm
        Y_aug[term_j, int_i] += -y_adm    # symmetric

    # The rest of off-diagonals already set from Ybus_full
    # Now we want Kron reduction so that we eliminate network buses and keep only internal nodes.
    keep_indices = list(range(0, M))  # retain only internal nodes
    try:
        Y_int_pre = kron_reduction(Y_aug, keep_indices)
        print("Kron reduction complete: Y_int_pre shape:", Y_int_pre.shape)
    except RuntimeError as e:
        print("Kron reduction failed:", e)
        Y_int_pre = None

    # If user requested a bus to be faulted (grounded), form the "faulted" physical matrix by removing that bus from network,
    # then redo Kron reduction (we must remove the bus row/col from network portion (indices M + (bus-1)) and adjust coupling matrices).
    Y_int_fault = None
    if faulted_bus_user is not None:
        b = int(faulted_bus_user)
        if not (1 <= b <= N):
            print("Requested faulted bus out of range, ignoring.")
        else:
            # remove bus b from the network part (indices M + (b-1))
            # build Y_aug_fault by deleting the corresponding row/col
            del_idx = M + (b-1)
            Y_aug_fault = np.delete(np.delete(Y_aug, del_idx, axis=0), del_idx, axis=1)
            # But note: we should also remove corresponding element from terminal bus connections; code above already does that by deleting
            try:
                Y_int_fault = kron_reduction(Y_aug_fault, keep_indices)
                print("Kron reduction on faulted case complete:", None if Y_int_fault is None else Y_int_fault.shape)
            except Exception as e:
                print("Kron fault reduction failed:", e)

    # If user requested a line trip (from_bus, to_bus), modify Ybus_full accordingly (remove off-diagonal and adjust diagonals)
    Y_int_post = None
    if line_trip is not None:
        a, b = line_trip
        if not (1 <= a <= N and 1 <= b <= N):
            print("Requested line trip endpoints out of range, ignoring.")
        else:
            # compute admittance of the line currently in Ybus_full: Y_ab = -Ybus_full[a-1, b-1]
            Y_line = -Ybus_full[a-1, b-1]
            print(f"Detected line admittance between {a} and {b} as {Y_line:.6g}")
            # Build modified physical matrix Y_mod where that off-diagonal is set to zero and both diagonals subtract Y_line
            Ybus_mod = Ybus_full.copy()
            # set off-diagonals
            Ybus_mod[a-1, b-1] = 0 + 0j
            Ybus_mod[b-1, a-1] = 0 + 0j
            # subtract diag contribution
            Ybus_mod[a-1, a-1] = Ybus_mod[a-1, a-1] - Y_line
            Ybus_mod[b-1, b-1] = Ybus_mod[b-1, b-1] - Y_line

            # rebuild augmented Y with same internal connections
            Y_aug_post = np.zeros((M+N, M+N), dtype=complex)
            Y_aug_post[M:M+N, M:M+N] = Ybus_mod
            # add internal coupling as before
            for idx, busnum in enumerate(gen_buses):
                terminal_idx = busnum - 1
                xd_p = float(xdp_list[idx])
                y_adm = -1j / xd_p
                int_i = idx
                term_j = M + terminal_idx
                Y_aug_post[int_i, int_i] += y_adm
                Y_aug_post[term_j, term_j] += y_adm
                Y_aug_post[int_i, term_j] += -y_adm
                Y_aug_post[term_j, int_i] += -y_adm
            try:
                Y_int_post = kron_reduction(Y_aug_post, keep_indices)
                print("Kron reduction post-line-trip complete.")
            except Exception as e:
                print("Kron post reduction failed:", e)

    # Save results
    save_dict = {
        "Ybus_full": Ybus_full,
        "bus_list": np.array(bus_list),
        "gen_buses": np.array(gen_buses),
        "Xd_primes": np.array(xdp_list),
        "Y_aug": Y_aug,
        "Y_int_pre": Y_int_pre
    }
    if Y_int_fault is not None:
        save_dict["Y_int_fault"] = Y_int_fault
    if Y_int_post is not None:
        save_dict["Y_int_post"] = Y_int_post

    sio.savemat(out_mat, save_dict)
    print("Saved matrices to:", out_mat)
    print("Done.")
