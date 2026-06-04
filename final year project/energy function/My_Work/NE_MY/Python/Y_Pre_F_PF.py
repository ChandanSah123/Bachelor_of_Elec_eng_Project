import os
import re
import numpy as np
import scipy.io as sio
import sys

# Attempt to import PSS/E modules
try:
    import psse3603 as psse_mod
    import psspy
    import dyntools
except ImportError:
    pass 

# ----------------- USER CONFIG -----------------
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\NE_MY"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_NE"
raw_file = os.path.join(work_dir, "IEEE39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
txt_file = os.path.join(result_dir, "Ybus_Export.txt")
mat_file = os.path.join(result_dir, "Y_all.mat")

# Fault & clearing specification (IEEE 39 Bus)
fault_bus = 29
trip_line_from = 29
trip_line_to = 26
trip_ckt_id = '1'

# ------------------------------------------------

def init_psse_and_read_case():
    psspy.psseinit(200)
    ierr = psspy.read(0, raw_file)
    if ierr != 0:
        raise RuntimeError(f"psspy.read returned error {ierr}")

def solve_and_convert():
    # 1. Solve Power Flow to get Voltages
    psspy.fnsl([0,0,0,1,1,0,99,0])
    
    # 2. Convert Generators Only (Norton Equivalent)
    psspy.cong(0)
    
    # 3. DO NOT CONVERT LOADS (conl)
    # We skip psspy.conl. We will calculate load admittance in Python.
    
    try:
        psspy.ordr(0)
        psspy.fact()
        psspy.tysl(0)
    except Exception:
        pass

def get_system_data_from_psse():
    """
    Extracts Generators, Loads, and Trip Line parameters directly from PSS/E memory.
    """
    print("Extracting System Data from PSS/E...")
    
    # --- 1. GENERATORS ---
    # Extract Gen Bus Numbers and X'd (Zsorce)
    ierr, mach_bus = psspy.amachint(-1, 4, 'NUMBER')
    ierr, mach_z   = psspy.amachcplx(-1, 4, 'ZSORCE')
    
    gen_buses = []
    xd_prime = []
    
    if mach_bus:
        for bus, z in zip(mach_bus[0], mach_z[0]):
            if abs(z) > 0:
                gen_buses.append(bus)
                xd_prime.append(z.imag) # X'd is the imaginary part
    
    print(f" > Found {len(gen_buses)} Generators.")

    # --- 2. LOADS (Calculate Y = S* / V^2) ---
    ierr, ld_bus = psspy.aloadint(-1, 4, ['NUMBER'])
    ierr, ld_mva = psspy.aloadcplx(-1, 4, ['MVAACT']) # Actual In-service Load
    ierr, bus_v  = psspy.abusreal(-1, 2, ['PU'])      # Bus Voltages
    
    # Create Voltage Map
    ierr, all_buses = psspy.abusint(-1, 2, ['NUMBER'])
    v_map = {b: v for b, v in zip(all_buses[0], bus_v[0])}

    load_adm = {}
    if ld_bus:
        for bus, s in zip(ld_bus[0], ld_mva[0]):
            v = v_map.get(bus, 1.0)
            if v < 0.001: v = 1.0 # Protect against bad voltage
            
            # Calculation: Y = conj(S_pu) / |V|^2
            s_pu = s / 100.0  # Convert MVA to p.u.
            y = np.conj(s_pu) / (v**2)
            
            # Add to dictionary (sum if multiple loads on one bus)
            load_adm[bus] = load_adm.get(bus, 0j) + y
            
    print(f" > Calculated Admittances for {len(load_adm)} Load Buses.")
    print(load_adm)

    # --- 3. TRIP LINE PARAMETERS ---
    # Find Y and B for the line 29-26
    y_trip = 0j
    b_trip = 0j
    
    ierr, br_int = psspy.abrnint(-1, 0, 0, 3, 2, ['FROMNUMBER', 'TONUMBER'])
    ierr, br_z   = psspy.abrncplx(-1, 0, 0, 3, 2, ['RX'])
    ierr, br_flt = psspy.abrnreal(-1, 0, 0, 3, 2, ['CHARG'])
    
    if br_int and br_z:
        # Handle case if CHARG is None (rare API bug)
        charges = br_flt[0] if (br_flt and br_flt[0]) else [0.0]*len(br_int[0])
        
        for f, t, z, ch in zip(br_int[0], br_int[1], br_z[0], charges):
            # Check match (Line 29-26 or 26-29)
            if (f == trip_line_from and t == trip_line_to) or (f == trip_line_to and t == trip_line_from):
                if abs(z) > 0:
                    y_trip = 1.0 / z
                    b_trip = 1j * ch / 2.0
                    print(f" > Found Trip Line {f}-{t}: Y={y_trip:.4f}, B_half={b_trip:.4f}")
                    break

    return gen_buses, xd_prime, load_adm, y_trip, b_trip

def export_ybus_to_text():
    ierr = psspy.output_y_matrix(0, 1, 0, 0, txt_file)
    if isinstance(ierr, tuple): ierr = ierr[0]
    if ierr != 0:
        raise RuntimeError(f"output_y_matrix error: {ierr}")

def parse_ybus_text_to_numpy():
    Yreal = {}
    Yimag = {}
    max_bus = 0
    float_re = re.compile(r"[-+]?\d*\.\d+|\d+")
    
    with open(txt_file, "r") as f:
        for line in f:
            parts = float_re.findall(line)
            if len(parts) >= 4:
                try:
                    i, j = int(parts[0]), int(parts[1])
                    real, imag = float(parts[2]), float(parts[3])
                    max_bus = max(max_bus, i, j)
                    Yreal[(i, j)] = real
                    Yimag[(i, j)] = imag
                except Exception:
                    continue
    
    # 39 Bus might have gaps or be sequential. We assume sequential 1..39 for matrix size.
    Y = np.zeros((max_bus, max_bus), dtype=np.complex128)
    for (i, j), rv in Yreal.items():
        iv = Yimag.get((i, j), 0.0)
        Y[i-1, j-1] = rv + 1j*iv
        Y[j-1, i-1] = rv + 1j*iv
    return Y

def build_faulted_physical(Y_full, fault_bus_num):
    """ Removes row/col of the fault bus """
    idx = fault_bus_num - 1
    Yf = np.delete(Y_full, idx, axis=0)
    Yf = np.delete(Yf, idx, axis=1)
    return Yf

def build_postfault_physical_exact(Y_full, a_bus, b_bus, y_line, b_charging):
    """ Removes line a-b using calculated parameters """
    idx_a = a_bus - 1
    idx_b = b_bus - 1
    Yp = Y_full.copy()
    
    # Subtract Y_series and B_charging from diagonals
    y_total = y_line + b_charging
    Yp[idx_a, idx_a] -= y_total
    Yp[idx_b, idx_b] -= y_total
    
    # Open the line (Zero off-diagonals)
    Yp[idx_a, idx_b] = 0.0 + 0.0j
    Yp[idx_b, idx_a] = 0.0 + 0.0j
    return Yp

def kron_reduce_internal(Y_physical, gen_buses_list, xd_prime_list, load_admittance=None):
    """
    Performs Kron Reduction. 
    """
    K = len(gen_buses_list)
    N = Y_physical.shape[0]

    # A. Calculate Generator Admittances (ybar)
    ybar = np.zeros((K, K), dtype=complex)
    for i, xd in enumerate(xd_prime_list):
        ybar[i, i] = 1.0 / (1j * xd)

    # B. Augment Physical Matrix (YD)
    YD = Y_physical.copy()
    
    # Detect if we are in Faulted Mode (Bus removed)
    # 39 Bus System: Pre=39x39, Fault=38x38
    max_bus_id = 39 # Hardcoded for IEEE 39
    is_faulted = (N < max_bus_id) 

    # Add Gen Admittances to Diagonals
    gen_indices_map = [] 
    for i, bus in enumerate(gen_buses_list):
        if is_faulted:
            if bus == fault_bus: continue 
            idx = bus - 1 if bus < fault_bus else bus - 2
        else:
            idx = bus - 1
            
        gen_indices_map.append((i, idx))
        YD[idx, idx] += ybar[i, i]

    # Add Load Admittances to Diagonals
    if load_admittance:
        for bus_num, yload in load_admittance.items():
            if is_faulted:
                if bus_num == fault_bus: continue
                idx = bus_num - 1 if bus_num < fault_bus else bus_num - 2
                YD[idx, idx] += yload
            else:
                idx = bus_num - 1
                YD[idx, idx] += yload

    # C. Build Partition Matrices
    YA = ybar.copy()
    YB = np.zeros((K, N), dtype=complex)
    
    for (k, phys_idx) in gen_indices_map:
        YB[k, phys_idx] = -ybar[k, k]

    YC = YB.T

    # D. Reduce
    try:
        YD_inv_YC = np.linalg.solve(YD, YC)
        Yint = YA - np.dot(YB, YD_inv_YC)
    except np.linalg.LinAlgError:
        print("Warning: Matrix singular.")
        Yint = np.zeros((K, K))
        
    return Yint

# === Main execution ===
if __name__ == "__main__":
    if not os.path.exists(raw_file):
        raise FileNotFoundError(f"raw file not found: {raw_file}")
    
    print("Initializing PSS/E...")
    init_psse_and_read_case()

    print("Solving...")
    solve_and_convert() 

    # --- GET DYNAMIC DATA (Generators, Loads, Line Params) ---
    gen_list, xd_list, load_adm, y_trip, b_trip = get_system_data_from_psse()

    print("Exporting Ybus (Network Only)...")
    export_ybus_to_text()
    
    print("Parsing Ybus...")
    Y_pre = parse_ybus_text_to_numpy()
    
    # -------------------------------------------------------------
    # 1. PRE-FAULT (With Calculated Loads)
    # -------------------------------------------------------------
    print("Calculating Pre-Fault Reduced Matrix...")
    Yint_pre = kron_reduce_internal(Y_pre, gen_list, xd_list, load_admittance=load_adm)

    # -------------------------------------------------------------
    # 2. FAULTED (Bus 29 Removed)
    # -------------------------------------------------------------
    print(f"Calculating Faulted Reduced Matrix (Bus {fault_bus} removed)...")
    Y_fault = build_faulted_physical(Y_pre, fault_bus)
    Yint_fault = kron_reduce_internal(Y_fault, gen_list, xd_list, load_admittance=load_adm)

    # -------------------------------------------------------------
    # 3. POST-FAULT (Line 29-26 Removed)
    # -------------------------------------------------------------
    print(f"Calculating Post-Fault Reduced Matrix (Line {trip_line_from}-{trip_line_to} removed)...")
    Y_post = build_postfault_physical_exact(Y_pre, trip_line_to, trip_line_from, y_trip, b_trip)
    Yint_post = kron_reduce_internal(Y_post, gen_list, xd_list, load_admittance=load_adm)

    # Verification Print
    #print("\n--- Python Result: Yint_Fault  ---")
    #print(Yint_fault)
    
    # Save
    print(f"\nSaving to {mat_file} ...")
    mdict = {
        "Y_pre": Y_pre,
        "Y_fault": Y_fault,
        "Y_post": Y_post,
        "Yint_pre": Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post": Yint_post,
        "gen_buses": gen_list  # Useful to know which row is which gen
    }
    sio.savemat(mat_file, mdict)
    print("SUCCESS: Matrices saved.")