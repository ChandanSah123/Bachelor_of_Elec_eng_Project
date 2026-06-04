%% 5. MOD Identification Algorithm (Iterative Lookup)
fprintf('\n==============================================\n');
fprintf('     STARTING MOD IDENTIFICATION LOOP\n');
fprintf('==============================================\n');

% Pre-calculate VPE (Potential Energy) as it does not depend on MOD
VPE = Calculate_PE(npts, g, Pi, C, D, th, ths);

% Prepare to store results
results_table = []; % Columns: [Num_Machines, CCT_TEF, Error, Vcr]
best_error = 9999;
best_MOD_group = [];
best_CCT_TEF = 0;

sorted_gen_indices = MOD_sort_data(:, 1); 
% LOOP: Try Top 1 machine, then Top 2 .....
for k = 1:num_gen
    fprintf('\n--- Testing MOD Candidate: Top %d Generator(s) ---\n', k);
    current_MOD_indices = sorted_gen_indices(1:k); 
    % 2. Calculate theta_u
    theta_u = Calculate_theta_u1(k, MOD_sort_data, num_gen, ths, H);
    % 3. Calculate CUEP for this specific MOD (Time Domain BCU)
    options = odeset('MaxStep',0.01,'InitialStep',0.01);
    [Tm, Ym] = ode45(@(t,y) Integrateth(t, y, Pi, C, D, H), [0 20], theta_u, options);
    theta_cuep_mod = Ym(end, :)';
    % 4. Calculate Critical Energy (Vcr) for this CUEP
    Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
    % 5. Calculate Corrected Kinetic Energy (Depends on MOD!)
    [~, KE_corr_candidate] = Calculate_KE(npts, g, H, Ws, w, current_MOD_indices);
    % 6. Total Energy & Margin
    V_total = VPE + KE_corr_candidate;
    delV_candidate = Vcr_candidate - V_total;
    % 7. Find CCT_TEF (Zero Crossing of Margin)
    cct_tef = NaN; % Default if no crossing found
    % Start search after fault inception
    idx_search = find(T >= 1.0, 1); 
    if isempty(idx_search), idx_search = 2; end
    for t = idx_search+1 : npts
        if delV_candidate(t) < 0 && delV_candidate(t-1) >= 0
            % Linear Interpolation
            y_prev = delV_candidate(t-1);
            y_curr = delV_candidate(t);
            fraction = y_prev / (y_prev - y_curr);
            cct_tef = T(t-1) + fraction * (T(t) - T(t-1));
            break; % Found the first crossing
        end
    end
    
    % 8. Error Calculation
    if isnan(cct_tef)
        fprintf('   -> System remained Stable (V < Vcr). No CCT found.\n');
        err = 9999;
    else
        err = abs(cct_tef - t_cct_absolute);
        fprintf('   -> CCT_TEF: %.4fs | CCT_TD: %.4fs | Error: %.4f\n', cct_tef, t_cct_absolute, err);
    end
    
    % Store in Table
    results_table = [results_table; k, cct_tef, err, Vcr_candidate];
    
    % Check if this is the best one
    if err < best_error
        best_error = err;
        best_MOD_group = current_MOD_indices;
        best_CCT_TEF = cct_tef;
    end
end

%% 7. BUILD/UPDATE OFFLINE DATABASE
% This block saves the "Learned" Information into a database file.

db_file = 'Offline_Database.mat';

% 1. Calculate the Kinetic Energy "Fingerprint" (Normalized KE) for this case
% We use the state at Fault Clearing Time (approx) or Max KE point
[~, idx_clear] = min(abs(T - (1.0 + 0.22))); % Example: Fault clears at 1.22s
w_online = w(idx_clear, :)';
KE_raw = 0.5 * M .* (w_online.^2);
KE_fingerprint = KE_raw / sum(KE_raw); % Normalized (Sum = 1)

% 2. Create the Entry Structure
new_entry.Fault_Location = 'Bus 9'; % You can automate this name
new_entry.MOD_Generators = best_MOD_group;
new_entry.CUEP_Angles    = theta_cuep_mod; % The CUEP we found using that MOD
new_entry.Critical_Energy= Vcr_candidate;  % The Vcr for that MOD
new_entry.KE_Signature   = KE_fingerprint; % The Search Key

% 3. Append to Database
if exist(db_file, 'file')
    load(db_file, 'TEF_Database');
    TEF_Database = [TEF_Database; new_entry];
else
    TEF_Database = new_entry;
end

save(db_file, 'TEF_Database');
fprintf('Entry added to Offline_Database.mat. Total Entries: %d\n', length(TEF_Database));