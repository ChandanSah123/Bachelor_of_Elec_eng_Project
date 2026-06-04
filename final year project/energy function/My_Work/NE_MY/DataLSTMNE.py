import psse3603
import psspy
import dyntools
import numpy as np
import pandas as pd
import os
import sys

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
WORK_DIR = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\NE"
RAW_FILE = os.path.join(WORK_DIR, "ieee39bus1.raw")
DYR_FILE = os.path.join(WORK_DIR, "ieee39buscls1.dyr")
OUT_FILE = os.path.join(WORK_DIR, "temp_master.out")
LOG_FILE = os.path.join(WORK_DIR, "psse_log.txt")

# Output Folder
DATASET_DIR = r"C:\Users\Acer\Desktop\LSTM_Master_Dataset"
os.makedirs(DATASET_DIR, exist_ok=True)

# --- SCENARIO PARAMETERS ---
# 1. Load Scalings (80%, 100%, 120%)
LOAD_SCALINGS = [80.0, 100.0, 120.0] 

# 2. Generator Outages (List of Gen Bus IDs to trip). 
# Use [0] for Base Case (No outage).
# Real Gen buses in IEEE39: 30-39.
GEN_OUTAGE_CASES = [0, 30, 38] # Example: Base Case, Trip Gen 30, Trip Gen 38. Add more if needed.

# 3. Fault Durations (Seconds)
# This sweep covers < CCT and > CCT for most buses.
FAULT_TIMES = [0.08, 0.12, 0.16, 0.22, 0.30] 

# 4. Generators to Record (Angles/Speeds)
GEN_BUSES = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39]

# PSS/E Initialization
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)
psspy.report_output(2, LOG_FILE, [0])

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

def get_all_buses():
    """Returns a list of all bus numbers in the case."""
    ierr, buses = psspy.abusint(-1, 2, "NUMBER")
    return buses[0]

def get_all_branches():
    """Returns a list of tuples (FromBus, ToBus) for all lines."""
    ierr, (from_b, to_b) = psspy.abrnint(-1, 0, 0, 2, 1, ["FROMNUMBER", "TONUMBER"])
    return list(zip(from_b, to_b))

def scale_load(percent):
    """Scales P and Q load by a percentage."""
    # 1. Subsystem selection (All buses)
    psspy.bsys(0, 0, [0.0, 0.0], 1, [1], 0, [], 0, [], 0, [])
    
    # 2. Scale Load (3 = constant MVA, constant current and constant admittance)
    # We scale Real (P) and Reactive (Q) power
    scale_factor = percent # e.g. 120.0
    psspy.scal_2(0, 0, 1, [0,0,0,0,0], [scale_factor, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

def trip_generator(gen_bus):
    """Disconnects a generator at the given bus."""
    if gen_bus == 0: return # Base Case
    # ID is usually '1' for single gen
    psspy.machine_data_2(gen_bus, r"""1""", [_i,_i,_i,_i,_i, 0], [_f]*17) # 0 = Out of service

# ==============================================================================
# 3. MAIN DATA GENERATION ENGINE
# ==============================================================================

print("--- Initializing Scenario Discovery ---")
# Load case once to get lists
psspy.read(0, RAW_FILE)
ALL_BUSES = get_all_buses() 
ALL_BRANCHES = get_all_branches()

# Reduce scope for testing? (Uncomment below to test small batch first)
# ALL_BUSES = [15, 29] 
# ALL_BRANCHES = [(16, 17)] 

print(f"Buses to Fault: {len(ALL_BUSES)}")
print(f"Lines to Fault: {len(ALL_BRANCHES)}")

sim_count = 0

# --- LOOP 1: LOAD LEVEL ---
for load_pct in LOAD_SCALINGS:
    
    # --- LOOP 2: TOPOLOGY (GEN OUTAGE) ---
    for out_gen in GEN_OUTAGE_CASES:
        
        # --- LOOP 3: FAULT DURATION ---
        for t_fault in FAULT_TIMES:
            
            # ---------------------------------------------
            # SCENARIO A: BUS FAULTS
            # ---------------------------------------------
            for f_bus in ALL_BUSES:
                sim_count += 1
                
                # 1. Reset Case
                psspy.read(0, RAW_FILE)
                psspy.dyre_new([1,1,1,1], DYR_FILE, "", "", "")
                
                # 2. Apply Conditions
                scale_load(load_pct)
                trip_generator(out_gen)
                
                # 3. Solve Power Flow (Important after changes!)
                ierr = psspy.fnsl([0,0,0,1,1,0,99,0])
                if ierr > 0:
                    print(f"Skipping Sim {sim_count} (Power Flow Diverged)")
                    continue

                # 4. Setup Dynamic Params
                psspy.dynamics_solution_param_2([_i]*8, [_f,_f, 0.001,_f,_f,_f,_f,_f])
                psspy.cong(0)
                psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
                psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
                psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
                
                # 5. Setup Channels (Angles & Speed)
                psspy.delete_all_plot_channels()
                for i, g_bus in enumerate(GEN_BUSES):
                    psspy.machine_array_channel([i+1, 1, g_bus], r"""1""", "") # Angle
                    psspy.machine_array_channel([i+11, 2, g_bus], r"""1""", "") # Speed

                # 6. Run Simulation (Bus Fault)
                psspy.strt(0, OUT_FILE)
                psspy.run(0, 1.0, 0, 1, 0)
                
                # Apply Bus Fault
                psspy.dist_bus_fault(f_bus, 1, 0.0, [0.0, -0.2E+10])
                
                psspy.run(0, 1.0 + t_fault, 0, 1, 0)
                psspy.dist_clear_fault(1)
                psspy.run(0, 5.0, 0, 1, 0)
                
                # 7. Extract & Save
                chnf = dyntools.CHNF(OUT_FILE)
                _, _, chandata = chnf.get_data()
                
                data_dict = {'Time': chandata['time']}
                for i, g_bus in enumerate(GEN_BUSES):
                    if (i+1) in chandata: data_dict[f'Gen{g_bus}_Ang'] = chandata[i+1]
                    if (i+11) in chandata: data_dict[f'Gen{g_bus}_Spd'] = chandata[i+11]
                
                df = pd.DataFrame(data_dict)
                
                # Filename: Ld[Load]_Out[Gen]_Bus[Bus]_T[Time].csv
                fname = f"Ld{int(load_pct)}_Out{out_gen}_Bus{f_bus}_T{t_fault}.csv"
                df.to_csv(os.path.join(DATASET_DIR, fname), index=False)
                
                # Progress Update
                if sim_count % 10 == 0: print(f"Saved {fname} (Total: {sim_count})")


            # ---------------------------------------------
            # SCENARIO B: LINE FAULTS
            # ---------------------------------------------
            for (from_b, to_b) in ALL_BRANCHES:
                sim_count += 1
                
                # 1. Reset Case
                psspy.read(0, RAW_FILE)
                psspy.dyre_new([1,1,1,1], DYR_FILE, "", "", "")
                
                # 2. Apply Conditions
                scale_load(load_pct)
                trip_generator(out_gen)
                
                # 3. Solve Power Flow
                ierr = psspy.fnsl([0,0,0,1,1,0,99,0])
                if ierr > 0: continue # Skip if bad case

                # 4. Dynamics Setup
                psspy.dynamics_solution_param_2([_i]*8, [_f,_f, 0.001,_f,_f,_f,_f,_f])
                psspy.cong(0)
                psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
                psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
                psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
                
                # 5. Channels
                psspy.delete_all_plot_channels()
                for i, g_bus in enumerate(GEN_BUSES):
                    psspy.machine_array_channel([i+1, 1, g_bus], r"""1""", "")
                    psspy.machine_array_channel([i+11, 2, g_bus], r"""1""", "")

                # 6. Run Simulation (Line Fault)
                psspy.strt(0, OUT_FILE)
                psspy.run(0, 1.0, 0, 1, 0)
                
                # Apply Line Fault (50% distance)
                # psspy.dist_branch_fault(from bus, to bus, ckt, percent)
                psspy.dist_branch_fault(from_b, to_b, r"""1""", 1, 0.0, [0.0, -0.2E+10])
                
                psspy.run(0, 1.0 + t_fault, 0, 1, 0)
                psspy.dist_clear_fault(1)
                psspy.run(0, 5.0, 0, 1, 0)
                
                # 7. Extract & Save
                chnf = dyntools.CHNF(OUT_FILE)
                _, _, chandata = chnf.get_data()
                
                data_dict = {'Time': chandata['time']}
                for i, g_bus in enumerate(GEN_BUSES):
                    if (i+1) in chandata: data_dict[f'Gen{g_bus}_Ang'] = chandata[i+1]
                    if (i+11) in chandata: data_dict[f'Gen{g_bus}_Spd'] = chandata[i+11]
                
                df = pd.DataFrame(data_dict)
                
                # Filename: Ld[Load]_Out[Gen]_Line[F-T]_T[Time].csv
                fname = f"Ld{int(load_pct)}_Out{out_gen}_Line{from_b}-{to_b}_T{t_fault}.csv"
                df.to_csv(os.path.join(DATASET_DIR, fname), index=False)

print(f"--- COMPLETE. Total Simulations: {sim_count} ---")