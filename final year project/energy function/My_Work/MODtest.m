% Original Data Setup
tic
g = 10;
ths = [-0.0199 0.2077 0.3073 0.2620 0.3868 0.2962 0.5595 0.2532 0.4642 -0.1625]; % Angles
H = [42 30.3 35.8 28.6 26 34.8 26.4 24.3 34.5 500]; % Weights/Heights
MOD_indices = 1:9; % Indices for Group 1 (The 'chain')

% 1. Grouping and Weighted Sums
H1 = sum(H(MOD_indices)); % Total weight of Group 1
H2 = H(g); % Total weight of Group 2 (Index 10)

th1 = sum(ths(MOD_indices) .* H(MOD_indices)) / H1; % Weighted average for Group 1
th2 = ths(g); % Weighted average for Group 2 (since H2 is just H(10) and th2 is just ths(10))

% 2. Calculating Correction Factors
dth = th1 - th2; % Difference
H_total = H1 + H2;

% Distribute the pi-based correction
dth1 = (pi - 2*dth) * (H2 / H_total); % Shift for Group 1
dth2 = (pi - 2*dth) * (H1 / H_total); % Shift for Group 2

% 3. Applying the Correction (Vectorized Update)
thu = ths; % Initialize the new angles vector

% Apply positive shift dth1 to indices 1-9
thu(MOD_indices) = ths(MOD_indices) + dth1;

% Apply negative shift dth2 to index 10
thu(g) = ths(g) - dth2;

% Output
disp('Original Angles (ths):')
disp(ths)
disp('Corrected Angles (thu):')
disp(thu)
toc