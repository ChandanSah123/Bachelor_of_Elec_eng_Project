%% Settings
t_fault = 1.0;      % Fault application time
t_cl = 0.168;       % Clearing time
t_end = 4.0;        % End time
del_t = 0.001; 
del_t_fault = 0.001; 
del_t_BCU = 0.01; 

num_bus = 39;
num_gen = 10;
plot_bus = 2; % Relative index for plotting (e.g., Bus 31)

% Directories (Update if needed)
% result_dir = '...'; 

%% Set up variables (Sorted by Bus 30 -> 39)
sys_case = 39;

% Load Data and Matrices
if ~exist('data', 'var'), load('data1.mat'); end
if ~exist('Yint_pre', 'var'), load('Y_all.mat'); end

% 1. Generator Parameters (Must be sorted Bus 30 -> 39)
% Bus: [30,   31,   32,   33,   34,   35,   36,   37,   38,   39]
H =    [42.0; 30.3; 35.8; 28.6; 26.0; 34.8; 26.4; 24.3; 34.5; 500.0];
xd_p = [0.031;0.0697;0.0531;0.0436;0.132;0.050;0.049;0.057;0.057;0.006];

% 2. Voltages (E) from Ma & Pai (Sorted 30->39)
E = [1.0562; 1.1782; 1.1507; 1.1004; 1.3594; 1.1806; 1.1288; 1.0721; 1.1382; 1.0215];

M_gen = H .* (1/(pi*60));
M_T = sum(M_gen);
D_gen = zeros(10,1);

idx_gen = [30 31 32 33 34 35 36 37 38 39];
idx_delta = 1:num_gen; 
idx_omega = num_gen+1:2*num_gen; 
num_var = 2*num_gen;

%% 3. Self-Consistent Initialization (THE FIX)
% We use the angles from PSS/E, but we CALCULATE Power to match Y-bus.
% Angles in PSS/E are Degrees -> Convert to Radians
angleE = data(1, 1:10)' * (pi/180); 

% Calculate Complex Voltage and Power Injection
V_pre = E .* exp(1j * angleE);
I_pre = Yint_pre * V_pre;
S_pre = V_pre .* conj(I_pre);

% FORCE Mechanical Power to equal Calculated Electrical Power
Pm_consistent = real(S_pre);
Pm = Pm_consistent; 

Pgen_pre = Pm; 
Pgen_post = Pm; 
Pgen = Pm;

disp('--- Initialization ---');
disp('Calculated Consistent Pm (p.u.):');
disp(Pm');

% Initial State
delpre = angleE;
del0 = sum(delpre .* M_gen) / M_T;
x_eq_pre = [delpre - del0; zeros(num_gen,1)]; 

Eeq_pre = E;
Eeq_post = E;

%% Kron Reduction Setup
Y_pre_kron = Yint_pre;
Y_fault_kron = Yint_fault;
Y_post_kron = Yint_post;

edge_kron = nchoosek(1:num_gen,2);
E_kron = zeros(size(edge_kron,1), num_gen);

% Pre-calculate vectors
gij_pre = zeros(size(edge_kron,1),1); bij_pre = zeros(size(edge_kron,1),1);
gij_flt = zeros(size(edge_kron,1),1); bij_flt = zeros(size(edge_kron,1),1);
gij_pst = zeros(size(edge_kron,1),1); bij_pst = zeros(size(edge_kron,1),1);

for i=1:size(E_kron,1)
    u = edge_kron(i,1); v = edge_kron(i,2);
    E_kron(i,u) = 1;  E_kron(i,v) = -1;
    
    prodE = Eeq_pre(u) * Eeq_pre(v);
    
    gij_pre(i) = real(Y_pre_kron(u,v))*prodE;  bij_pre(i) = imag(Y_pre_kron(u,v))*prodE;
    gij_flt(i) = real(Y_fault_kron(u,v))*prodE; bij_flt(i) = imag(Y_fault_kron(u,v))*prodE;
    gij_pst(i) = real(Y_post_kron(u,v))*prodE;  bij_pst(i) = imag(Y_post_kron(u,v))*prodE;
end

P_pre_kron   = Pgen_pre  - real(diag(Y_pre_kron))  .* Eeq_pre.^2;
P_fault_kron = Pgen_post - real(diag(Y_fault_kron)).* Eeq_post.^2;
P_post_kron  = Pgen_post - real(diag(Y_post_kron)) .* Eeq_post.^2;

%% Model Functions
f_pre_kron= @(x) [x(idx_omega)-1/M_T*M_gen'*x(idx_omega); diag(1./M_gen)*(-D_gen.*x(idx_omega)-E_kron'*diag(bij_pre)*sin(E_kron*x(idx_delta))-abs(E_kron)'*diag(gij_pre)*cos(E_kron*x(idx_delta))+P_pre_kron)-1/M_T*(sum(P_pre_kron)-2*ones(size(E_kron'))*diag(gij_pre)*cos(E_kron*x(idx_delta)))];
J_pre_kron=@(x) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_pre)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_pre)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_pre)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];

f_fault_kron= @(x) [x(idx_omega)-1/M_T*M_gen'*x(idx_omega); diag(1./M_gen)*(-D_gen.*x(idx_omega)-E_kron'*diag(bij_flt)*sin(E_kron*x(idx_delta))-abs(E_kron)'*diag(gij_flt)*cos(E_kron*x(idx_delta))+P_fault_kron)-1/M_T*(sum(P_pre_kron)-2*ones(size(E_kron'))*diag(gij_flt)*cos(E_kron*x(idx_delta)))];
J_fault_kron=@(x) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_flt)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_flt)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_flt)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];

f_post_kron= @(x) [x(idx_omega)-1/M_T*M_gen'*x(idx_omega); diag(1./M_gen)*(-D_gen.*x(idx_omega)-E_kron'*diag(bij_pst)*sin(E_kron*x(idx_delta))-abs(E_kron)'*diag(gij_pst)*cos(E_kron*x(idx_delta))+P_post_kron)-1/M_T*(sum(P_post_kron)-2*ones(size(E_kron'))*diag(gij_pst)*cos(E_kron*x(idx_delta)))];
J_post_kron=@(x) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_pst)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_pst)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_pst)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];

f_post_theta= @(x) -E_kron'*diag(bij_pst)*sin(E_kron*x)-abs(E_kron)'*diag(gij_pst)*cos(E_kron*x)+P_post_kron-(M_gen/M_T).*(sum(P_post_kron)-2*ones(size(E_kron'))*diag(gij_pst)*cos(E_kron*x));
J_post_theta=@(x) -E_kron'*diag(bij_pst)*diag(cos(E_kron*x))*E_kron+abs(E_kron)'*diag(gij_pst)*diag(sin(E_kron*x))*E_kron-2*diag(M_gen)/M_T*ones(size(E_kron'))*diag(gij_pst)*diag(sin(E_kron*x))*E_kron;

%% ---- Solve post-fault SEP using ROBUST NR ----
disp('Calculating Post-Fault SEP...');
theta01 = x_eq_pre(idx_delta); 
[theta_sep, converged] = RobustNR(f_post_theta, J_post_theta, theta01);

if ~converged
    error('Post-fault SEP did not converge even with Robust Solver.');
end
disp('SEP Converged!');
x_eq_post = [theta_sep; zeros(num_gen,1)];

%% Direct Simulation (TDS)
tic; 
num_sim_idx = ceil((t_end-t_cl)/del_t) + ceil(t_cl/del_t_fault) + 100;
x_sim = zeros(num_var, num_sim_idx); 
t_sim = zeros(1, num_sim_idx);
x_sim(:,1) = x_eq_pre; t_sim(1)=0; t_idx=1;

while t_sim(t_idx) <= t_end
    if t_sim(t_idx) < t_fault
        f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_pre_kron(x_sim(:,t_idx))+f_pre_kron(x));
        J=@(x) -eye(num_var)+del_t/2*J_pre_kron(x);
        dt=del_t;
    elseif t_sim(t_idx) < t_fault+t_cl
        f=@(x) -x+x_sim(:,t_idx)+del_t_fault/2*(f_fault_kron(x_sim(:,t_idx))+f_fault_kron(x));
        J=@(x) -eye(num_var)+del_t_fault/2*J_fault_kron(x);
        dt=del_t_fault;
    else
        f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_post_kron(x_sim(:,t_idx))+f_post_kron(x));
        J=@(x) -eye(num_var)+del_t/2*J_post_kron(x);
        dt=del_t;
    end
    
    [x_next, ok] = RobustNR(f, J, x_sim(:,t_idx));
    if ~ok, disp(['TDS Failed at t=' num2str(t_sim(t_idx))]); break; end
    
    x_sim(:,t_idx+1) = x_next;
    t_sim(t_idx+1) = t_sim(t_idx) + dt;
    t_idx = t_idx + 1;
end
t_sim=t_sim(1:t_idx); x_sim=x_sim(:,1:t_idx);
tcomp_TDS=toc;
result_TDS=(max(x_sim(idx_delta,t_idx-1))-min(x_sim(idx_delta,t_idx-1)))<2*pi;

%% PEBS Method
tic; x_fault=x_eq_pre;
max_PEBS_sim=30000; % Increased limit
V_ke=zeros(1,max_PEBS_sim); V_p=zeros(1,max_PEBS_sim); V_d=zeros(1,max_PEBS_sim); V_total=zeros(1,max_PEBS_sim);
PEBS=zeros(1,max_PEBS_sim);

PEBS(1) = f_post_theta(x_fault(idx_delta,1))'*(x_fault(idx_delta,1)-x_eq_post(idx_delta));

V_pe0 = -P_post_kron'*(x_eq_pre(idx_delta)-x_eq_post(idx_delta)) ...
        -bij_pst'*(cos(E_kron*x_eq_pre(idx_delta))-cos(E_kron*x_eq_post(idx_delta))) ...
        +0.5*(gij_pst.*(cos(E_kron*x_eq_pre(idx_delta))+cos(E_kron*x_eq_post(idx_delta))))' * ...
        (abs(E_kron)*x_eq_pre(idx_delta)-abs(E_kron)*x_eq_post(idx_delta));

PEBS_idx=0; t_idx=1;
idx_tcl=round(t_cl/del_t_fault);

while t_idx < max_PEBS_sim && (PEBS_idx==0 || t_idx <= idx_tcl + 500)
    f=@(x) -x+x_fault(:,t_idx)+del_t_fault/2*(f_fault_kron(x_fault(:,t_idx))+f_fault_kron(x));
    J=@(x) -eye(num_var)+del_t_fault/2*J_fault_kron(x);
    
    [x_next, ok] = RobustNR(f,J,x_fault(:,t_idx));
    if ~ok, break; end
    x_fault(:,t_idx+1) = x_next;
    
    % Energy Calculation
    curr_theta = x_fault(idx_delta,t_idx+1);
    
    V_ke(t_idx+1) = 0.5*sum(M_gen.*x_fault(idx_omega,t_idx+1).^2);
    V_p(t_idx+1) = -P_post_kron'*(curr_theta-x_eq_post(idx_delta)) ...
                   -bij_pst'*(cos(E_kron*curr_theta)-cos(E_kron*x_eq_post(idx_delta)));
    V_d(t_idx+1) = V_d(t_idx)+0.5*(gij_pst.*(cos(E_kron*curr_theta)+cos(E_kron*x_fault(idx_delta,t_idx))))' * ...
                   (abs(E_kron)*curr_theta - abs(E_kron)*x_fault(idx_delta,t_idx));
    V_total(t_idx+1) = V_ke(t_idx+1) + V_p(t_idx+1) + V_d(t_idx+1);
    
    PEBS(t_idx+1) = f_post_theta(curr_theta)'*(curr_theta - x_eq_post(idx_delta));
    
    if PEBS(t_idx) < 0 && PEBS(t_idx+1) >= 0 && PEBS_idx == 0
        PEBS_idx = t_idx; 
    end
    t_idx=t_idx+1;
end
t_fault_sim=(0:t_idx-1)*del_t_fault; 

% Handle Missing Crossing
if PEBS_idx == 0
    warning('PEBS Crossing not found. Simulation may be too short or system too stable.');
    PEBS_idx = t_idx; % Prevents Array Index Error
end

Vcr_PEBS = V_p(PEBS_idx) + V_d(PEBS_idx) - V_pe0; 
idx_tcr_PEBS_arr = find(V_total(1:t_idx) < Vcr_PEBS); 
if isempty(idx_tcr_PEBS_arr)
    tcr_PEBS = 0;
else
    tcr_PEBS = del_t_fault * idx_tcr_PEBS_arr(end);
end

tcomp_PEBS=toc;

%% Display Results
disp(['PEBS CCT: ' num2str(tcr_PEBS)]);

% ==========================================================
% HELPER FUNCTION: ROBUST NEWTON RAPHSON (Replaces old NR)
% ==========================================================
function [x, converged] = RobustNR(f_func, J_func, x0)
    x = x0;
    converged = false;
    max_iter = 50;
    damping = 0.8;
    
    for k = 1:max_iter
        F = f_func(x);
        err = norm(F, inf);
        if err < 1e-4
            converged = true;
            return;
        end
        
        J = J_func(x);
        
        % Robust Solver for Singular/Ill-Conditioned Jacobian
        if rcond(J) < 1e-12
            dx = -pinv(J) * F;
        else
            dx = -J \ F;
        end
        
        x = x + damping * dx;
    end
end