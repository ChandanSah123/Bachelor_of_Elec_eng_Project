import psse3603
import psspy
import dyntools
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io
import numpy as np
import os
import re

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC programs"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")
dyr_file = os.path.join(work_dir, "ieee9bus.dyr")
out_file = os.path.join(work_dir, "IEEE9.out")
excel_file = os.path.join(work_dir, "IEEE9_Results.xlsx")
mat_file = os.path.join(work_dir, "IEEE9_Results.mat")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation for IEEE 9-Bus ---")

# Simulation Parameters
fault_bus = 7
from_bus=7
to_bus=5
line_id='1'
fault_time = 0
clear_time = 1
end_time = 4

# Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# Solver Settings
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# --- DEFINE CHANNELS (UPDATED) ---
# NOTE: Assumes Generators are at Bus 1, 2, and 3 with ID '1'

#Channels

psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,6,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])

# --- CONVERT NETWORK ---
print("Converting Network...")
psspy.cong(0) 

# Load Conversion (3 Steps)
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])

# Run Simulation
print("Initializing State (STRT)...")
ierr_strt = psspy.strt(0, out_file)

if ierr_strt > 0:
    print(f"ERROR: Simulation failed to start. Error code: {ierr_strt}")
else:
    psspy.run(0, fault_time, 0, 1, 0)

    print(f"Applying Fault at Bus {fault_bus}...")
    psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])

    psspy.run(0, clear_time, 0, 1, 0)
    print("Clearing Fault by TRIPPING LINE 7-5...")
    psspy.dist_branch_trip(from_bus, to_bus, line_id)   # Use correct line ID

    print("Clearing Fault...")
    psspy.dist_clear_fault(1)

    psspy.run(0, end_time, 0, 1, 0)
    print("Simulation Complete.")

# ==============================================================================
# 3. DATA EXTRACTION
# ==============================================================================

if not os.path.exists(out_file):
    print(f"CRITICAL ERROR: Output file not found at {out_file}")
else:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
    time_data = chandata['time']

# ==============================================================================
# 5. EXPORT
# ==============================================================================
print("Preparing Export...")
t_data = chandata['time']
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])

# Prepare a single DataFrame for Excel export
final_df = pd.DataFrame({'Time': t_data})  # Start with a single Time column

for key in sorted_keys:
    if key in chandata:
        y_data = chandata[key]
        raw_name = str(chanid[key]).strip()
        final_df[raw_name] = y_data  # Add only the data column

# Save the DataFrame to Excel
try:
    final_df.to_excel(excel_file, index=False)
    print(f"Excel saved to {excel_file}")
except:
    print("Could not save Excel (Check if file is open).")

# MATLAB Export
mat_export = {}
data1_export = []  # To store data for data1.mat

for key in sorted_keys:
    if key in chandata:
        y_data = chandata[key]
        raw_name = str(chanid[key]).strip()

        # MATLAB
        clean_name = re.sub(r'[^a-zA-Z0-9]', '_', raw_name)
        clean_name = re.sub(r'_+', '_', clean_name).strip('_')
        if clean_name and clean_name[0].isdigit():
             clean_name = "Ch_" + clean_name
        min_len = min(len(t_data), len(y_data))
        matrix_nx2 = np.column_stack((t_data[:min_len], y_data[:min_len]))
        mat_export[clean_name] = matrix_nx2

        # Collect data for data1.mat (exclude time column)
        data1_export.append(y_data[:min_len])

# Save ieee9_result.mat
try:
    if mat_export:
        scipy.io.savemat(mat_file, mat_export)
        print(f"MAT file saved to {mat_file}")
except Exception as e:
    print(f"Error saving MAT: {e}")

# Save data1.mat
try:
    if data1_export:
        data1_matrix = np.column_stack(data1_export)  # Combine all data columns
        data1_file = os.path.join(work_dir, "data1.mat")
        scipy.io.savemat(data1_file, {"data": data1_matrix})
        print(f"data1.mat file saved to {data1_file}")
except Exception as e:
    print(f"Error saving data1.mat: {e}")

print("--- Script Finished ---")