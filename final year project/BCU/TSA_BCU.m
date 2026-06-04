%% TSA_BCU.m  - refactored for readability (keeps original algorithm and results)
% Main script implementing:
%  - Kron-reduction-based TEF (pre/fault/post)
%  - Time-domain simulation (TDS)
%  - PEBS computation
%  - BCU (Boundary Controlling UEP) method
%
% Original author: Dongchan Lee (August 2017)
% Refactored for clarity and readability; algorithm unchanged.

clear; close all;

%% ------------------ USER SETTINGS ------------------
sys_case = 9;           % IEEE system (calls dyn9.m)
fault_line = 4;         % index of faulted line in Line.con
fault_frto_bus = 1;     % use from(1) or to(2) endpoint of the line for fault
t_cl = 0.15;            % clearing time (s)

% Simulation time settings
t_end = 20;             % simulation end time (s)
del_t = 0.05;           % main simulation timestep (s)
del_t_fault = 0.001;    % time step during fault interval (s)
del_t_BCU = 0.01;       % timestep for BCU iteration (s)
t_fault = 0;            % fault application time (s)

% Plotting selection
plot_bus = 2;

%% ------------------ LOAD SYSTEM DATA ------------------
run(['dyn' int2str(sys_case)])    % loads Bus.con, Line.con, Syn.con, etc.

% Indices and sizes
slack_bus = SW.con(1);
num_bus = size(Bus.con,1);
num_line = size(Line.con,1);
num_gen = size(Syn.con,1);
num_load = num_bus - num_gen;

% Index helper variables
eye_bus = eye(num_bus);
idx_gen = Syn.con(:,1);
idx_load = setdiff(1:num_bus, idx_gen);
idx_delta = 1:num_gen;
idx_omega = num_gen + 1 : 2*num_gen;
num_var = 2 * num_gen;

% Generator parameters
% M_gen is moment of inertia (inertia-like) used in the code: Syn.con(:,18)/(2*pi*freq)
M_gen = Syn.con(:,18) ./ (2*pi*Syn.con(:,4));   % consistent with original code
M_T = sum(M_gen);                               % total inertia-like sum
Syn.con(:,19) = 2 * ones(size(Syn.con(:,19))); % set damping to 2 (keeps original script behavior)
D_gen = Syn.con(:,19) ./ (2*pi*Syn.con(:,4));

% Voltage setpoints and d-axis transient reactance
v_gen = [SW.con(:,4); PV.con(:,5)];  % V for slack and PV buses
xd_p = Syn.con(:,9);                 % Xd' values

% Line / network construction
line_frto = Line.con(:,1:2);                             % from-to mapping
fault_bus = line_frto(fault_line, fault_frto_bus);       % faulted bus index
Z_line = Line.con(:,8) + 1i * Line.con(:,9);             % series impedance (complex)
E = eye_bus(line_frto(:,1), :) - eye_bus(line_frto(:,2), :); % incidence matrix for lines
Y = E' * diag((Z_line) .^ -1) * E;                       % network admittance

% Generator and load injections
Pgen = zeros(num_bus,1);
Pgen(PV.con(:,1)) = PV.con(:,4);                         % generation P at PV buses
Sload = zeros(num_bus,1);
Sload(PQ.con(:,1)) = (PQ.con(:,4) + 1i * PQ.con(:,5));   % load complex S at PQ buses

%% ------------------ Convert static loads to constant impedance form (ZIP -> Z) -----------
% Solve powerflow-like equation to get V_eq, then compute y_load = conj(S)/V^2
x_eq = NR_ss(Y, Pgen - Sload, idx_load, v_gen, slack_bus);   % solve S_inj = Pgen - Sload
V_eq = x_eq(num_bus+1:end) .* (cos(x_eq(1:num_bus)) + 1i * sin(x_eq(1:num_bus))); 
I_eq = Y * V_eq;
S_inj = V_eq .* conj(I_eq);
% Compute equivalent admittances of loads at operating point (to use as Z loads)
y_load = conj(Sload) ./ (V_eq .^ 2);

% Build network Y matrices for pre / fault / post situations
YN_pre = E' * diag(Z_line .^ -1) * E + diag(y_load);
YN_fault = YN_pre;
Zf_line = Z_line; Zf_line(fault_line) = inf;   % fault: set line impedance to infinite (open) for Y_post construction
YN_post = E' * diag(Zf_line .^ -1) * E + diag(y_load);

%% ------------------ Build augmented network including machine stator impedances ----------
% augmented dimension = num_bus + num_gen
Y_pre = zeros(num_bus + num_gen);
% Place stator block and coupling between internal machine nodes and network nodes
% Partitioned as [internal Machines ; network nodes]
Y_pre([1:num_gen, num_gen + idx_gen'], [1:num_gen, num_gen + idx_gen']) = ...
    [ diag( (1i * xd_p) .^ -1 ),  diag( -(1i * xd_p) .^ -1 );
      diag( -(1i * xd_p) .^ -1 ), diag(  (1i * xd_p) .^ -1 ) ];

% Add network YN_pre to appropriate block
Y_pre(num_gen+1:end, num_gen+1:end) = Y_pre(num_gen+1:end, num_gen+1:end) + YN_pre;

% build Y_fault and Y_post analogously (copy then modify)
Y_fault = Y_pre; Y_post = Y_pre;
Y_fault(num_gen+1:end, num_gen+1:end) = Y_fault(num_gen+1:end, num_gen+1:end) + YN_fault;
% remove the faulted network node (corresponds to a 3-phase fault on that bus)
Y_fault(num_gen + fault_bus, :) = []; Y_fault(:, num_gen + fault_bus) = [];
% for post-fault
Y_post(num_gen+1:end, num_gen+1:end) = Y_post(num_gen+1:end, num_gen+1:end) + YN_post;

%% ------------------ Compute pre/post equilibrium and internal EMFs --------------------
% Pre-contingency equilibrium (reduced network YN_pre used for load modelling)
x_eq_pre = NR_ss(YN_pre, Pgen, idx_load, v_gen, slack_bus);
V_eq_pre = x_eq_pre(num_bus+1:end) .* (cos(x_eq_pre(1:num_bus)) + 1i * sin(x_eq_pre(1:num_bus)));
I_eq_pre = YN_pre * V_eq_pre;
Pgen_pre = real( V_eq_pre(idx_gen) .* conj(I_eq_pre(idx_gen)) );
% Internal EMF magnitude and angle (E' = V + j Xd' * I)
Eeq_pre = abs( V_eq_pre(idx_gen) + 1i * xd_p .* I_eq_pre(idx_gen) );
delta_eq_pre = angle( V_eq_pre(idx_gen) + 1i * xd_p .* I_eq_pre(idx_gen) );
x_eq_pre = [ delta_eq_pre - (M_gen' * delta_eq_pre) / M_T; zeros(num_gen,1) ];   % COI-referenced initial state

% Post-contingency equilibrium (network after clearing = YN_post)
x_eq_post = NR_ss(YN_post, Pgen, idx_load, v_gen, slack_bus);
V_eq_post = x_eq_post(num_bus+1:end) .* (cos(x_eq_post(1:num_bus)) + 1i * sin(x_eq_post(1:num_bus)));
I_eq_post = YN_post * V_eq_post;
Pgen_post = real( V_eq_post(idx_gen) .* conj(I_eq_post(idx_gen)) );
Eeq_post = abs( V_eq_post(idx_gen) + 1i * xd_p .* I_eq_post(idx_gen) );
delta_eq_post = angle( V_eq_post(idx_gen) + 1i * xd_p .* I_eq_post(idx_gen) );
x_eq_post = [ delta_eq_post - (M_gen' * delta_eq_post) / M_T; zeros(num_gen,1) ];

%% ------------------ Kron Reduction to generator internal nodes (pre, fault, post) ------
% Kron reduce Y_pre / Y_fault / Y_post to the generator internal nodes (1:num_gen)
Y_pre_kron = kron_reduce(Y_pre, num_gen);
Y_fault_kron = kron_reduce(Y_fault, num_gen);
Y_post_kron = kron_reduce(Y_post, num_gen);

% Build pairwise edge incidence (for E matrix style used in subsequent energy computations)
edge_kron = nchoosek(1:num_gen, 2);
E_kron = zeros(size(edge_kron,1), num_gen);
for i = 1:size(edge_kron,1)
    E_kron(i, edge_kron(i,1)) = 1;
    E_kron(i, edge_kron(i,2)) = -1;
end

% Compute gij and bij (scaled by internal EMFs) for each pair (pre / fault / post)
gij_pre_kron = zeros(size(E_kron,1),1);
gij_fault_kron = zeros(size(E_kron,1),1);
gij_post_kron = zeros(size(E_kron,1),1);
bij_pre_kron = zeros(size(E_kron,1),1);
bij_fault_kron = zeros(size(E_kron,1),1);
bij_post_kron = zeros(size(E_kron,1),1);

for i = 1:size(E_kron,1)
    a = edge_kron(i,1); b = edge_kron(i,2);
    gij_pre_kron(i)   = real(Y_pre_kron(a,b)) * Eeq_pre(a) * Eeq_pre(b);
    gij_fault_kron(i) = real(Y_fault_kron(a,b)) * Eeq_post(a) * Eeq_post(b);
    gij_post_kron(i)  = real(Y_post_kron(a,b)) * Eeq_post(a) * Eeq_post(b);
    bij_pre_kron(i)   = imag(Y_pre_kron(a,b)) * Eeq_pre(a) * Eeq_pre(b);
    bij_fault_kron(i) = imag(Y_fault_kron(a,b)) * Eeq_post(a) * Eeq_post(b);
    bij_post_kron(i)  = imag(Y_post_kron(a,b)) * Eeq_post(a) * Eeq_post(b);
end

% Compute P_i terms used inside potential energy (P_i = Pm_i - E_i^2 * G_ii)
P_pre_kron   = Pgen_pre - real(diag(Y_pre_kron)) .* (Eeq_pre .^ 2);
P_fault_kron = Pgen_post - real(diag(Y_fault_kron)) .* (Eeq_post .^ 2);
P_post_kron  = Pgen_post - real(diag(Y_post_kron)) .* (Eeq_post .^ 2);

%% ------------------ Define model RHS and Jacobian (pre, fault, post) ---------------
% To improve readability we create small helper functions that compute f(x) and J(x)
% f_pre_kron(x), J_pre_kron(x) etc. will be used by the time-stepping solver and NR().

% For clarity: x = [delta; omega] with COI reference applied consistent with x_eq_pre/post
f_pre_kron = @(x) model_rhs(x, E_kron, bij_pre_kron, gij_pre_kron, P_pre_kron, M_gen, M_T, D_gen);
J_pre_kron = @(x) model_jacobian(x, E_kron, bij_pre_kron, gij_pre_kron, P_pre_kron, M_gen, M_T, D_gen);

f_fault_kron = @(x) model_rhs(x, E_kron, bij_fault_kron, gij_fault_kron, P_fault_kron, M_gen, M_T, D_gen);
J_fault_kron = @(x) model_jacobian(x, E_kron, bij_fault_kron, gij_fault_kron, P_fault_kron, M_gen, M_T, D_gen);

f_post_kron = @(x) model_rhs(x, E_kron, bij_post_kron, gij_post_kron, P_post_kron, M_gen, M_T, D_gen);
J_post_kron = @(x) model_jacobian(x, E_kron, bij_post_kron, gij_post_kron, P_post_kron, M_gen, M_T, D_gen);

% f_post_theta and J_post_theta used in BCU for theta-only dynamics
f_post_theta = @(theta) theta_rhs(theta, E_kron, bij_post_kron, gij_post_kron, P_post_kron, M_gen, M_T);
J_post_theta = @(theta) theta_jacobian(theta, E_kron, bij_post_kron, gij_post_kron, P_post_kron, M_gen, M_T);

%% ------------------ Time Domain Simulation (Direct Integration using implicit trapezoid + NR) -----
tic;
num_sim_idx = ceil((t_end - t_cl)/del_t) + ceil(t_cl / del_t_fault);
x_sim = zeros(num_var, num_sim_idx);
t_sim = zeros(1, num_sim_idx);

x_sim(:,1) = x_eq_pre;  % start at pre-contingency equilibrium
t_sim(1) = 0; t_idx = 1;

while t_sim(t_idx) <= t_end
    % Choose appropriate model depending on time window (pre-fault, fault-on, post-fault)
    if t_sim(t_idx) < t_fault
        f_handle = @(x) -x + x_sim(:,t_idx) + del_t/2 * ( f_pre_kron(x_sim(:,t_idx)) + f_pre_kron(x) );
        J_handle = @(x) -eye(num_var) + del_t/2 * J_pre_kron(x);
        t_sim(t_idx+1) = t_sim(t_idx) + del_t;
    elseif t_sim(t_idx) < t_fault + t_cl
        f_handle = @(x) -x + x_sim(:,t_idx) + del_t_fault/2 * ( f_fault_kron(x_sim(:,t_idx)) + f_fault_kron(x) );
        J_handle = @(x) -eye(num_var) + del_t_fault/2 * J_fault_kron(x);
        t_sim(t_idx+1) = t_sim(t_idx) + del_t_fault;
    else
        f_handle = @(x) -x + x_sim(:,t_idx) + del_t/2 * ( f_post_kron(x_sim(:,t_idx)) + f_post_kron(x) );
        J_handle = @(x) -eye(num_var) + del_t/2 * J_post_kron(x);
        t_sim(t_idx+1) = t_sim(t_idx) + del_t;
    end

    % Solve implicit step by Newton (NR wrapper)
    x_sim(:, t_idx+1) = NR(f_handle, J_handle, x_sim(:, t_idx));
    t_idx = t_idx + 1;
end

% Trim arrays to actual simulated length
t_sim = t_sim(1:t_idx);
x_sim = x_sim(:, 1:t_idx);
tcomp_TDS = toc;

% result_TDS true if angles bounded (not runaway)
result_TDS = ( max( x_sim(idx_delta, end) ) - min( x_sim(idx_delta, end) ) ) < 2*pi;

%% ------------------ PEBS Computation (Potential Energy Boundary Surface) -----------------
tic;
x_fault = x_eq_pre;    % initial state for fault-on simulation
max_PEBS_sim = 10000;

V_ke = zeros(1, max_PEBS_sim);
V_p  = zeros(1, max_PEBS_sim);
V_d  = zeros(1, max_PEBS_sim);

% initial energies
V_ke(1) = 0.5 * sum( M_gen .* x_fault(idx_omega,1) .^ 2 );
V_p(1)  = - P_post_kron' * ( x_fault(idx_delta,1) - x_eq_post(idx_delta) ) - bij_post_kron' * ( cos(E_kron * x_fault(idx_delta,1)) - cos(E_kron * x_eq_post(idx_delta)) );
V_d(1)  = 0;

% baseline potential energy at t=0 (used in comparisons)
V_pe0 = - P_post_kron' * ( x_eq_pre(idx_delta) - x_eq_post(idx_delta) ) ...
        - bij_post_kron' * ( cos(E_kron * x_eq_pre(idx_delta)) - cos(E_kron * x_eq_post(idx_delta)) ) ...
        + 0.5 * ( (gij_post_kron .* ( cos(E_kron * x_eq_pre(idx_delta)) + cos(E_kron * x_eq_post(idx_delta)) ))' * ( abs(E_kron) * x_eq_pre(idx_delta) - abs(E_kron) * x_eq_post(idx_delta) ) );

V_total = V_ke + V_p + V_d;
PEBS = zeros(1, max_PEBS_sim);
PEBS(1) = f_post_theta( x_fault(idx_delta,1) )' * ( x_fault(idx_delta,1) - x_eq_post(idx_delta) );

PEBS_idx = 0; t_idx = 1;
idx_tcl = round( (t_cl + t_fault) / del_t_fault );
idx_tf = max(round(t_fault / del_t_fault), 1);

% Iterate fault-on trajectory (implicit trapezoid)
while (t_idx < max_PEBS_sim) && ( (PEBS_idx == 0) || (t_idx <= idx_tcl) )
    f_handle = @(x) -x + x_fault(:,t_idx) + del_t_fault/2 * ( f_fault_kron(x_fault(:,t_idx)) + f_fault_kron(x) );
    J_handle = @(x) -eye(num_var) + del_t_fault/2 * J_fault_kron(x);
    x_fault(:, t_idx+1) = NR(f_handle, J_handle, x_fault(:, t_idx));

    % update energies
    V_ke(t_idx+1) = 0.5 * sum( M_gen .* x_fault(idx_omega, t_idx+1) .^ 2 );
    V_p(t_idx+1)  = - P_post_kron' * ( x_fault(idx_delta, t_idx+1) - x_eq_post(idx_delta) ) - bij_post_kron' * ( cos(E_kron * x_fault(idx_delta, t_idx+1)) - cos(E_kron * x_eq_post(idx_delta)) );
    V_d(t_idx+1)  = V_d(t_idx) + 0.5 * ( (gij_post_kron .* ( cos(E_kron * x_fault(idx_delta, t_idx+1)) + cos(E_kron * x_fault(idx_delta, t_idx)) ))' * ( abs(E_kron) * x_fault(idx_delta, t_idx+1) - abs(E_kron) * x_fault(idx_delta, t_idx) ) );

    V_total(t_idx+1) = V_ke(t_idx+1) + V_p(t_idx+1) + V_d(t_idx+1);
    PEBS(t_idx+1) = f_post_theta( x_fault(idx_delta, t_idx+1) )' * ( x_fault(idx_delta, t_idx+1) - x_eq_post(idx_delta) );

    % detect PEBS crossing (zero crossing)
    if (PEBS(t_idx) < 0) && (PEBS(t_idx+1) >= 0) && (PEBS_idx == 0)
        PEBS_idx = t_idx;
    end
    t_idx = t_idx + 1;
end

t_fault_sim = (0:t_idx-1) * del_t_fault;
% trim arrays
V_ke = V_ke(1:t_idx); V_p = V_p(1:t_idx); V_d = V_d(1:t_idx); V_total = V_total(1:t_idx); PEBS = PEBS(1:t_idx);

% compute critical energy and times
Vcr_PEBS = V_p(PEBS_idx) + V_d(PEBS_idx) - V_pe0;
idx_tcr_PEBS = find(V_total < Vcr_PEBS); idx_tcr_PEBS = idx_tcr_PEBS(end);
tcr_PEBS = del_t_fault * idx_tcr_PEBS;
result_PEBS = V_total(idx_tcl) < Vcr_PEBS;
tcomp_PEBS = toc;

%% ------------------ BCU (Boundary Controlling UEP) method --------------------------
tic;
t_idx = 1; BCU_idx = 0;
max_BCU_sim = round(5 / del_t_BCU);

theta_sim = zeros(num_gen, max_BCU_sim);
theta_abs = zeros(1, max_BCU_sim);

% Initialization for BCU: start from the PEBS-detected state
theta_sim(:,1) = x_fault(idx_delta, PEBS_idx);
theta_abs(1) = sum( abs( f_post_theta( theta_sim(:,1) ) ) );

while (t_idx < max_BCU_sim) && (PEBS_idx ~= 0) && (BCU_idx == 0)
    % Trapezoidal step for theta-only dynamics
    f_handle = @(x) -x + theta_sim(:,t_idx) + del_t_BCU/2 * ( f_post_theta(theta_sim(:,t_idx)) + f_post_theta(x) );
    J_handle = @(x) -eye(num_gen) + del_t_BCU/2 * J_post_theta(x);

    theta_sim(:, t_idx+1) = NR(f_handle, J_handle, theta_sim(:, t_idx));
    theta_abs(t_idx+1) = sum( abs( f_post_theta( theta_sim(:, t_idx+1) ) ) );

    t_idx = t_idx + 1;

    % heuristic: detect peak in |f| that indicates uep
    if ( theta_abs(t_idx) > theta_abs(t_idx-1) ) && ( theta_abs(max(1, t_idx-2)) > theta_abs(t_idx-1) )
        BCU_idx = t_idx - 2;
    end
end

t_bcu_sim = (0:t_idx-1) * del_t_BCU;
theta_sim = theta_sim(:, 1:t_idx); theta_abs = theta_abs(:, 1:t_idx);

% Unstable equilibrium found by BCU
theta_u = theta_sim(:, BCU_idx);
x_eq_unstable = [ theta_u; zeros(num_gen,1) ];

% Compute Vcr using BCU formula from the paper (keeps original algebra)
Vcr_BCU = - P_post_kron' * ( theta_u - x_eq_post(idx_delta) ) ...
          - bij_post_kron' * ( cos(E_kron * theta_u) - cos(E_kron * x_eq_post(idx_delta)) ) ...
          + ( ( abs(E_kron) * (theta_u - x_eq_post(idx_delta)) ) ./ ( E_kron * (theta_u - x_eq_post(idx_delta)) ) )' * diag(gij_post_kron) * ( sin(E_kron * theta_u) - sin(E_kron * x_eq_post(idx_delta)) ) ...
          - V_pe0;

% find BCU-based tcr
idx_tcr_BCU = find(V_total < Vcr_BCU); idx_tcr_BCU = idx_tcr_BCU(end);
tcr_BCU = del_t_fault * idx_tcr_BCU;
result_BCU = V_total(idx_tcl) < Vcr_BCU;
tcomp_BCU = toc; tcomp_BCU = tcomp_BCU + tcomp_PEBS;

%% ------------------ DISPLAY / PLOTS (same as original) ----------------------------
% (plotting and display preserved)
figure;
subplot(2,1,1); hold on; grid on; box on;
plot(t_sim, x_sim(idx_delta,:));
xlim([0 t_end]); set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('\delta');
subplot(2,1,2); hold on; grid on; box on;
plot(t_sim, x_sim(idx_omega,:));
xlim([0 t_end]); set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('\omega');

figure;
subplot(3,1,1); hold on; grid on; box on;
plot(t_fault_sim, V_ke, 'r--');
plot(t_fault_sim, V_p + V_d, 'b--');
plot(t_fault_sim, V_total);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('Energy function');
legend('V_{ke}','V_{pe}','V_{ke}+V_{pe}');

subplot(3,1,2); hold on; grid on; box on;
plot(t_fault_sim, PEBS);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('PEBS function');

subplot(3,1,3); hold on; grid on; box on;
plot(t_bcu_sim, theta_abs);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('iteration'); ylabel('|f|');

% Phase portrait and energy contour (preserved plotting approach)
figure;
plot_axis = [plot_bus plot_bus + num_gen];
plot_rng = [-pi pi -10 10];
x_plot_pre = x_eq_pre;
x_plot_post = x_eq_post;
[x1_plot_ij, x2_plot_ij] = meshgrid(linspace(plot_rng(1), plot_rng(2), 30), linspace(plot_rng(3), plot_rng(4), 30));
V_plot = zeros(size(x1_plot_ij));
dx1_plot_pre = zeros(size(x1_plot_ij)); dx2_plot_pre = dx1_plot_pre;
dx1_plot_post = dx1_plot_pre; dx2_plot_post = dx1_plot_pre;

for i = 1:size(x1_plot_ij,1)
    for j = 1:size(x1_plot_ij,2)
        x_plot_pre(plot_axis) = [ x1_plot_ij(i,j); x2_plot_ij(i,j) ];
        x_plot_post(plot_axis)= [ x1_plot_ij(i,j); x2_plot_ij(i,j) ];
        dx_plot_pre = f_fault_kron(x_plot_pre);
        dx_plot_post = f_post_kron(x_plot_post);
        dx1_plot_pre(i,j) = dx_plot_pre(plot_axis(1));
        dx2_plot_pre(i,j) = dx_plot_pre(plot_axis(2));
        dx1_plot_post(i,j) = dx_plot_post(plot_axis(1));
        dx2_plot_post(i,j) = dx_plot_post(plot_axis(2));
        num_prevent_zero = E_kron * ( x_plot_post(idx_delta) - x_eq_post(idx_delta) );
        num_prevent_zero(num_prevent_zero == 0) = inf;
        V_plot(i,j) = 0.5 * sum( M_gen .* x_plot_post(idx_omega) .^ 2 ) ...
                      - P_post_kron' * ( x_plot_post(idx_delta) - x_eq_post(idx_delta) ) ...
                      - bij_post_kron' * ( cos(E_kron * x_plot_post(idx_delta)) - cos(E_kron * x_eq_post(idx_delta)) ) ...
                      + ( ( abs(E_kron) * (x_plot_post(idx_delta) - x_eq_post(idx_delta)) ) ./ ( num_prevent_zero ) )' * diag(gij_post_kron) * ( sin(E_kron * x_plot_post(idx_delta)) - sin(E_kron * x_eq_post(idx_delta)) );
    end
end

subplot(1,2,1); hold on; grid on; box on;
plot( x_sim(plot_axis(1), idx_tf:idx_tcl), x_sim(plot_axis(2), idx_tf:idx_tcl), 'LineWidth', 3 );
streamslice(x1_plot_ij, x2_plot_ij, dx1_plot_pre, dx2_plot_pre);
axis(plot_rng);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); title('Fault-on Dynamics');

subplot(1,2,2); hold on; grid on; box on;
scatter(theta_u(plot_axis(1)), 0, 100, 'r', 'filled');
scatter(x_fault(plot_axis(1), PEBS_idx), 0, 100, 'b', 'filled');
plot(x_sim(plot_axis(1), :), x_sim(plot_axis(2), :), 'b', 'LineWidth', 2);
streamslice(x1_plot_ij, x2_plot_ij, dx1_plot_post, dx2_plot_post);
axis(plot_rng);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); title('Post-contingency Dynamics'); legend('u.e.p. (BCU)','u.e.p. (PEBS)');

% Energy surface plot
figure; hold on; grid on; box on;
scatter(theta_u(plot_axis(1)), 0, 100, 'r', 'filled');
scatter(x_fault(plot_axis(1), PEBS_idx), 0, 100, 'b', 'filled');
streamslice(x1_plot_ij, x2_plot_ij, dx1_plot_post, dx2_plot_post);
contour(x1_plot_ij, x2_plot_ij, V_plot, 'LineWidth', 2);
mesh(x1_plot_ij, x2_plot_ij, V_plot);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); zlabel('V'); legend('u.e.p. (BCU)','u.e.p. (PEBS)');

%% ------------------ Final display summary -------------------
disp(['System: ' num2str(sys_case) ' bus system with fault on line ' num2str(fault_line) ' and bus ' num2str(line_frto(fault_line, fault_frto_bus))]);
disp(['Result: (TDS) ' repmat('un',1,~result_TDS) 'stable  (BCU) ' repmat('un',1,~result_BCU) 'stable  (PEBS) ' repmat('un',1,~result_PEBS) 'stable']);
disp(['V critical: (@t_cl) ' num2str(V_total(idx_tcl),'%.3f') ' (BCU) ' num2str(Vcr_BCU,'%.3f') ' (PEBS) ' num2str(Vcr_PEBS,'%.3f')]);
disp(['Critical Clearing Time: (t_cl) ' num2str(t_cl,'%.3f') ' (BCU) ' num2str(tcr_BCU,'%.3f') ' (PEBS) ' num2str(tcr_PEBS,'%.3f')]);
disp(['Computation Time: (TDS) ' num2str(tcomp_TDS,'%.2f') ' (BCU) ' num2str(tcomp_BCU,'%.2f') ' (PEBS) ' num2str(tcomp_PEBS,'%.2f')]);

%% ------------------ Helper functions (local to this file) -------------------------
% Note: these helpers preserve the exact algebra of the original script.

function Yk = kron_reduce(Yfull, nGen)
% Kron reduction of Yfull to the first nGen nodes (internal machines)
% Yfull is partitioned as:
% [ Y11  Y12 ;
%   Y21  Y22 ]
Y11 = Yfull(1:nGen, 1:nGen);
Y12 = Yfull(1:nGen, nGen+1:end);
Y21 = Yfull(nGen+1:end, 1:nGen);
Y22 = Yfull(nGen+1:end, nGen+1:end);
% Reduced Y
Yk = Y11 - Y12 * (Y22 \ Y21);
end

function rhs = model_rhs(x, E_kron, bij_vec, gij_vec, P_vec, Mgen, MT, Dgen)
% model_rhs: builds the RHS of the full (theta, omega) system
% x: [theta; omega]
n = length(Mgen);
m = size(E_kron,1);

theta = x(1:n);
omega = x(n+1:2*n);

% first eq: angle derivative (omega with COI adjustment)
% Note: (Mgen' * omega) is scalar, subtract from each element
angle_dot = omega - (1/MT) * (Mgen' * omega);

% intermediate vectors (m x 1)
sin_vec = sin(E_kron * theta);   % m x 1
cos_vec = cos(E_kron * theta);   % m x 1

% terms (n x 1)
term_sin = E_kron' * ( diag(bij_vec) * sin_vec );             % n x 1
term_cos = abs(E_kron)' * ( diag(gij_vec) * cos_vec );       % n x 1

% COI correction: ones(n,m) * ( diag(gij_vec) * cos_vec ) -> n x 1
COI_corr = ones(n, m) * ( diag(gij_vec) * cos_vec );         % n x 1

% Build acceleration (n x 1)
acc = (1 ./ Mgen) .* ( - Dgen .* omega - term_sin - term_cos + P_vec ) ...
      - (1/MT) * ( sum(P_vec) - 2 * COI_corr );

rhs = [ angle_dot; acc ];
end

function J = model_jacobian(x, E_kron, bij_vec, gij_vec, P_vec, Mgen, MT, Dgen)
% model_jacobian: builds Jacobian corresponding to model_rhs
n = length(Mgen);
m = size(E_kron,1);
theta = x(1:n);

% Zero/Identity blocks
J11 = zeros(n);
J12 = eye(n) - ones(n,1) * (Mgen' / MT);

% build trig vectors
cos_vec = cos(E_kron * theta);   % m x 1
sin_vec = sin(E_kron * theta);   % m x 1

% bottom-left: derivative of acc wrt theta (n x n)
% A_theta = - E_kron' * diag(bij_vec) * diag(cos_vec) * E_kron
A_theta = - E_kron' * ( diag(bij_vec) * diag(cos_vec) ) * E_kron;

% B_theta = + abs(E_kron)' * diag(gij_vec) * diag(sin_vec) * E_kron
B_theta = abs(E_kron)' * ( diag(gij_vec) * diag(sin_vec) ) * E_kron;

% COI term derivative: (2/MT) * ones(n,m) * ( diag(gij_vec) * diag(sin_vec) ) * E_kron
COI_term = (2 / MT) * ( ones(n, m) * ( diag(gij_vec) * diag(sin_vec) ) ) * E_kron;

J21 = diag(1 ./ Mgen) * ( A_theta + B_theta ) - COI_term;

% bottom-right: derivative wrt omega
J22 = - diag(Dgen ./ Mgen);

% assemble Jacobian
J = [ J11, J12; J21, J22 ];
end

function val = theta_rhs(theta, E_kron, bij_vec, gij_vec, P_vec, Mgen, MT)
% RHS for theta-only BCU iteration (used by f_post_theta)
n = length(Mgen);
m = size(E_kron,1);

sin_vec = sin(E_kron * theta);    % m x 1
cos_vec = cos(E_kron * theta);    % m x 1

term1 = - E_kron' * ( diag(bij_vec) * sin_vec );          % n x 1
term2 = - abs(E_kron)' * ( diag(gij_vec) * cos_vec );    % n x 1
COI_corr = ones(n, m) * ( diag(gij_vec) * cos_vec );     % n x 1

val = term1 + term2 + P_vec - (Mgen / MT) .* ( sum(P_vec) - 2 * COI_corr );
end

function Jt = theta_jacobian(theta, E_kron, bij_vec, gij_vec, P_vec, Mgen, MT)
% Jacobian of theta-only RHS
n = length(Mgen);
m = size(E_kron,1);

cos_vec = cos(E_kron * theta);   % m x 1
sin_vec = sin(E_kron * theta);   % m x 1

% -E_kron' * diag(bij) * diag(cos) * E_kron
A = - E_kron' * ( diag(bij_vec) * diag(cos_vec) ) * E_kron;

% + abs(E_kron)' * diag(gij) * diag(sin) * E_kron
B = abs(E_kron)' * ( diag(gij_vec) * diag(sin_vec) ) * E_kron;

% COI related term: 2 * diag(Mgen) / MT * ones(n,m) * ( diag(gij) * diag(sin) ) * E_kron
COI_mat = 2 * diag(Mgen) / MT * ( ones(n, m) * ( diag(gij_vec) * diag(sin_vec) ) ) * E_kron;

Jt = A + B - COI_mat;
end
