import psse3603
import psspy
import dyntools
import os

# ==============================================================================
# 1. SETUP
# ==============================================================================
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_NE\Python Script"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF Implementation\39bus"
raw_file = os.path.join(work_dir, "IEEE39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
out_file = os.path.join(result_dir, "IEEE9_CCT_Search.out")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION FUNCTION
# ==============================================================================
def run_simulation(clearing_duration):
    """
    Returns: (is_stable, max_spread_found)
    """
    # --- Match Parameters Exactly to Verification Script ---
    fault_bus = 4
    from_bus = 4
    to_bus = 14
    line_id = '1 '
    
    t_fault_start = 1.0
    t_clear = t_fault_start + clearing_duration
    t_end = 4.0  # CHANGED to 4.0 to match your verification script

    # --- Clean Previous Run ---
    if os.path.exists(out_file):
        try:
            os.remove(out_file)
        except:
            pass

    # --- Load & Convert ---
    psspy.read(0, raw_file)
    psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

    # Use 0.01s Time Step
    psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
    psspy.fnsl([0,0,0,1,1,0,99,0])

    psspy.cong(0) 
    psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.fact()
    psspy.tysl(0)

    # --- Channels ---
    psspy.delete_all_plot_channels()
    psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0]) # Angles

    # --- Run ---
    psspy.strt(0, out_file)
    psspy.run(0, t_fault_start, 0, 1, 0)
    
    # Fault
    psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])
    psspy.run(0, t_clear, 0, 1, 0)
    
    # Clear
    psspy.dist_branch_trip(from_bus, to_bus, line_id)
    psspy.dist_clear_fault(1)
    psspy.run(0, t_end, 0, 1, 0)

    # --- Analysis ---
    try:
        chnf_obj = dyntools.CHNF(out_file)
        _, chanid, chandata = chnf_obj.get_data()
    except:
        return False, 999.0
    
    if 'time' not in chandata or len(chandata['time']) == 0:
        return False, 999.0

    angle_keys = [k for k in chanid.keys() if isinstance(k, int)]
    num_steps = len(chandata['time'])
    
    max_spread_over_run = 0.0
    is_stable = True
    
    for i in range(num_steps):
        angles = [chandata[k][i] for k in angle_keys]
        if angles:
            spread = max(angles) - min(angles)
            if spread > max_spread_over_run:
                max_spread_over_run = spread
            
            # STABILITY THRESHOLD
            # If your verification plot shows spread > 180 but stable, 
            # increase this to 360.
            if spread > 180.0: 
                is_stable = False
                # We don't break immediately so we can report the TRUE max spread
    
    return is_stable, max_spread_over_run

# ==============================================================================
# 3. MAIN CCT SEARCH LOOP
# ==============================================================================
print("\n========================================")
print("   Starting CCT Search (Verbose)")
print("========================================\n")

psspy.report_output(2, "", [])

t_min = 0.1  # Start closer to expected range to save time
t_max = 1
precision = 0.001 

print(f"Searching between {t_min}s and {t_max}s...\n")

while (t_max - t_min) > precision:
    t_test = (t_max + t_min) / 2
    
    is_stable, spread = run_simulation(t_test)
    
    status = "STABLE  " if is_stable else "UNSTABLE"
    print(f"t_cl = {t_test:.4f} s -> {status} (Max Spread: {spread:.2f} deg)")
    
    if is_stable:
        t_min = t_test 
    else:
        t_max = t_test 
        
print("\n========================================")
print(f"  CRITICAL CLEARING TIME (CCT): {t_min:.4f} s")
print("========================================")