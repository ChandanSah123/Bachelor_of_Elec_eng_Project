import psse3603
import psspy
import dyntools
import scipy.io
import numpy as np
import os

# ==============================================================================
# 1. SETUP
# ==============================================================================
# Update these paths to your IEEE 39 folder
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_Pre\39_bus\Python script"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_Pre\39_bus\Matlab script"

raw_file = os.path.join(work_dir, "IEEE39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
out_file = os.path.join(result_dir, "IEEE39.out")
mat_file = os.path.join(result_dir, "data1.mat") 

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation for IEEE 39-Bus ---")

# --- Simulation Settings ---
# We use a clearing time (0.2s) slightly HIGHER than CCT (0.185s)
# to ensure the system goes unstable for MOD identification.
fault_bus = 4
from_bus = 4
to_bus = 14
line_id = '1 '

fault_time = 1.0
t_cl = 0.20          # Clearing duration (must be > CCT for offline analysis)
clear_time = fault_time + t_cl
end_time = 4.0

# --- Load Case ---
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# --- Solve Power Flow ---
# Critical: Use 0.001s step size for accuracy
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# ==============================================================================
# 3. CHANNEL SETUP (10 GENERATORS)
# ==============================================================================
# The order here determines the column order in data1.mat
# IEEE 39 has 10 Generators. The [-1, -1, -1] selects ALL of them.
psspy.delete_all_plot_channels()

# 1. ANGLE (Code 1) -> Becomes Channels 1-10
psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0])

# 2. SPEED (Code 7) -> Becomes Channels 11-20
psspy.chsb(0, 1, [-1, -1, -1, 1, 7, 0])

# 3. MECHANICAL POWER (Code 6) -> Becomes Channels 21-30
psspy.chsb(0, 1, [-1, -1, -1, 1, 6, 0])

# 4. ELECTRICAL POWER (Code 2) -> Becomes Channels 31-40
psspy.chsb(0, 1, [-1, -1, -1, 1, 2, 0])

# 5. REACTIVE POWER (Code 4) -> Becomes Channels 41-50
psspy.chsb(0, 1, [-1, -1, -1, 1, 4, 0])

# ==============================================================================
# 4. CONVERT & RUN
# ==============================================================================
print("Converting Network...")
psspy.cong(0) 
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.fact()
psspy.tysl(0)

# Run
ierr = psspy.strt(0, out_file)
if ierr != 0:
    print(f"Error starting simulation: {ierr}")
    exit()

psspy.run(0, fault_time, 0, 1, 0)

print(f"Applying Fault at Bus {fault_bus}...")
psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, clear_time, 0, 1, 0)

print(f"Tripping Line {from_bus}-{to_bus}...")
psspy.dist_branch_trip(from_bus, to_bus, line_id)
psspy.dist_clear_fault(1)

psspy.run(0, end_time, 0, 1, 0)

# ==============================================================================
# 5. EXTRACTION & EXPORT
# ==============================================================================
print("Extracting Data...")
try:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
except Exception as e:
    print(f"Error extracting data: {e}")
    exit()

t_data = chandata['time']

# Filter only integer keys (channel IDs) and sort them
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])

data_columns = []

# Loop to build the matrix
for key in sorted_keys:
    y_data = chandata[key]
    min_len = min(len(t_data), len(y_data))
    data_columns.append(y_data[:min_len])

# Stack columns horizontally
if data_columns:
    final_matrix = np.column_stack(data_columns)
    
    # Save to .mat
    scipy.io.savemat(mat_file, {"data": final_matrix})
    
    print(f"SUCCESS: {mat_file} saved.")
    print(f"Matrix Shape: {final_matrix.shape}")
    print("-" * 30)
    print("DATA STRUCTURE MAPPING:")
    print("Columns 1-10 : Rotor Angle (delta)")
    print("Columns 11-20: Speed (omega)")
    print("Columns 21-30: Mechanical Power (Pm)")
    print("Columns 31-40: Electrical Power (Pe)")
    print("Columns 41-50: Reactive Power (Q)")
    print("-" * 30)
else:
    print("Error: No data extracted.")