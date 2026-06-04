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
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\NE_MY"
result_dir=r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\NE_MY\Result_Directory"
raw_file = os.path.join(work_dir, "ieee39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
out_file = os.path.join(result_dir, "NE.out")
excel_file = os.path.join(result_dir,"Simulation_Results.xlsx")
mat_file =os.path.join(result_dir, "Simulation_Results.mat")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation ---")

# Simulation Parameters
tf = [0.87, 0.328, 0.306, 0.311, 0.298, 0.29, 0.35, 0.342, 1, 0.293, 0.309, 1, 0.318, 0.314, 0.276, 0.203, 0.242, 0.293, 0.245, 0.284, 0.238, 0.21, 0.197, 0.256, 0.27, 0.213, 0.292, 0.215, 0.183, 1, 0.323, 0.294, 0.272, 0.275, 0.262, 0.189, 0.292, 0.176]
b = 22
fault_time = 1
clear_time = 1.3
end_time = 10

# Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# Solver Settings
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# Define Channels (What to record)
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,3,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])

# Convert Network
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])

# Run Simulation
psspy.strt(0, out_file)
psspy.run(0, fault_time, 0, 1, 0)

print(f"Applying Fault at Bus {b}...")
psspy.dist_bus_fault(b, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, clear_time, 0, 1, 0)

print("Clearing Fault...")
psspy.dist_clear_fault(1)

psspy.run(0, end_time, 0, 1, 0)
print("Simulation Complete.")

# ==============================================================================
# 3. DATA EXTRACTION
# ==============================================================================
chnf_obj = dyntools.CHNF(out_file)
short_title, chanid, chandata = chnf_obj.get_data()
time_data = chandata['time']
# ==============================================================================
# 5. EXPORT (Time attached to EVERY column)
# ==============================================================================
print("Preparing Export...")

# --- Common Setup ---
t_data = chandata['time']

# FIX IS HERE: Filter keys so we only sort Integers (Channel Numbers)
# This ignores 'time' string or other metadata keys that cause the TypeError
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])

# Storage for Excel (List of DataFrames)
excel_dfs = []

# Storage for MATLAB (Dictionary)
mat_export = {}
data1_export = [] 

print("Processing channels...")
for key in sorted_keys:
    # Double check key is in chandata
    if key in chandata:
        # Get Data
        y_data = chandata[key]
        
        # -------------------------------------------------
        # A. Prepare for EXCEL (DataFrames)
        # -------------------------------------------------
        raw_name = str(chanid[key]).strip()
        
        # Create mini-table: [Time, Value]
        mini_df = pd.DataFrame({
            'Time': t_data,
            raw_name: y_data
        })
        excel_dfs.append(mini_df)

        # -------------------------------------------------
        # B. Prepare for MATLAB (2-Column Matrices)
        # -------------------------------------------------
        # 1. Sanitize Name
        clean_name = re.sub(r'[^a-zA-Z0-9]', '_', raw_name)
        clean_name = re.sub(r'_+', '_', clean_name).strip('_')
        
        if clean_name[0].isdigit():
            clean_name = "Ch_" + clean_name
        min_len = min(len(t_data), len(y_data))

        # 2. Create an Nx2 Matrix [Time, Value]
        matrix_nx2 = np.column_stack((t_data, y_data))
        
        # 3. Add to dictionary
        mat_export[clean_name] = matrix_nx2

           # Collect data for data1.mat (exclude time column)
        data1_export.append(y_data[:min_len])

# --- Save Excel ---
if excel_dfs:
    final_df = pd.concat(excel_dfs, axis=1)
    try:
        final_df.to_excel(excel_file, index=False)
        print(f"SUCCESS: Excel saved to {excel_file}")
    except PermissionError:
        print("ERROR: Excel file is OPEN. Please close it.")

# --- Save MATLAB ---
try:
    if mat_export:
        scipy.io.savemat(mat_file, mat_export)
        print(f"SUCCESS: MAT file saved to {mat_file}")
    else:
        print("Warning: No data found to export to MATLAB.")
except Exception as e:
    print(f"Error saving MAT file: {e}")

# Save data1.mat
try:
    if data1_export:
        data1_matrix = np.column_stack(data1_export)  # Combine all data columns
        data1_file = os.path.join(result_dir, "data1.mat")
        scipy.io.savemat(data1_file, {"data": data1_matrix})
        print(f"data1.mat file saved to {data1_file}")
except Exception as e:
    print(f"Error saving data1.mat: {e}")

print("--- Script Finished ---")