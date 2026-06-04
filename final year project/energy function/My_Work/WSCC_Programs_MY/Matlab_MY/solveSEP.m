function [x_coi, x_raw, exitflag, lambda, resid] = solveSEP(Yint_post, Pm, E, H, x0)
% solveSEP  Solve for the Stable Equilibrium Point (SEP) of a multi-machine
%           post-fault system using constrained optimization.
%
% USAGE
%   [x_coi, x_raw, exitflag, lambda, resid] = solveSEP(Yint_post, Pm, E, H, x0)
%
% INPUTS
%   Yint_post : (n x n) complex internal bus admittance matrix (post-fault)
%   Pm        : (1 x n) mechanical powers (pu)
%   E         : (1 x n) internal voltage magnitudes (behind reactance)
%   H         : (1 x n) inertia constants (or M_i weights)
%   x0        : (1 x n) initial guess for rotor angles (rad)
%
% OUTPUTS
%   x_coi  : (1 x n) rotor angles in COI reference (rad)
%   x_raw  : (1 x n) raw rotor angles returned by the optimizer (rad)
%   exitflag: fmincon exit flag
%   lambda : Lagrange multipliers structure from fmincon
%   resid  : infinity norm of equality constraint residuals (should be ~0)
%
% NOTES
%   - This solver uses fmincon to find x such that ceq(x)=0 where ceq enforces
%     Pm = Pe (with COI correction). The objective is a dummy constant.
%   - Angles are constrained to [-pi, pi] for numerical stability.
%   - The final result is shifted to COI: x_coi = x_raw - delta_coi.
%
% Example (use your known parameters):
%   H = [23.64 6.4 3.01];
%   E = [1.0728 1.0775 1.0609];
%   Pm = [0.7195 1.63 0.85];
%   Yint_post = Y1; % supply your post-fault internal Y matrix
%   x0 = data(100,1:3); % or any initial guess
%   [x_coi, x_raw, exitflag, lambda, resid] = solveSEP(Yint_post, Pm, E, H, x0);

%% Input checks
n = numel(Pm);
assert(all(size(Yint_post) == [n n]), 'Yint_post must be n x n');
assert(numel(E) == n && numel(H) == n && numel(x0) == n, 'Pm, E, H, x0 must be length n');

%% Precompute C and D matrices (coupling coefficients)
C = zeros(n,n); % multiplies sin(delta_i - delta_j)
D = zeros(n,n); % multiplies cos(delta_i - delta_j)
for i=1:n
    for j=1:n
        C(i,j) = E(i)*E(j) * imag(Yint_post(i,j));
        D(i,j) = E(i)*E(j) * real(Yint_post(i,j));
    end
end

%% Optimization setup
% Bounds: keep angles within [-pi, pi]
lb = -pi * ones(1,n);
ub =  pi * ones(1,n);

% Empty linear constraints
A = []; b = []; Aeq = []; beq = [];

% fmincon options
opts = optimset('Algorithm','interior-point', ...
                'Display','off', ...
                'TolFun',1e-10, ...
                'TolX',1e-10, ...
                'MaxIter',1000);

% Dummy objective function: constant (we're solving via constraints)
obj = @(x) 0;

% Nonlinear constraint handle (calls nested function nlc_sep below)
nonlcon = @(x) nlc_sep(x, Pm, E, C, D, H, Yint_post);

%% Solve using fmincon
try
    [x_raw, ~, exitflag, output, lambda] = fmincon(obj, x0, A,b,Aeq,beq, lb, ub, nonlcon, opts);
catch ME
    % If fmincon not available or error, rethrow with friendly message
    rethrow(ME);
end

%% Normalize angles into [-pi,pi] (wrap)
x_raw = wrapToPi(x_raw);

%% Shift to COI frame
delta_coi = sum(x_raw .* H) / sum(H);
x_coi = x_raw - delta_coi;

%% Constraint residuals (sanity)
[~, ceq] = nlc_sep(x_raw, Pm, E, C, D, H, Yint_post);
resid = norm(ceq, Inf);

%% Print basic diagnostics if not converged well
if exitflag <= 0 || resid > 1e-6
    warning('solveSEP:exitOrResid','fmincon exitflag=%d, residual=%g', exitflag, resid);
end

end

%% ---------------- Nested helper function ----------------
function [c, ceq] = nlc_sep(x, Pm, E, C, D, H, Y1)
    % nlc_sep Compute nonlinear equality constraints for SEP in COI frame.
    %
    % Returns:
    %   c   : [] (no inequality)
    %   ceq : vector of length n containing Pm_i - Pe_i(x) with COI correction
    %
    % This mirrors the algebraic equilibrium used in TEF/COI formulations.

    c = []; % no inequality constraints
    nloc = numel(Pm);

    % Ensure angle wrapping to [-pi,pi] for numerical stability
    x = wrapToPi(x);

    % Compute Pcoi = sum_{i<j} 2*D(i,j)*cos(xi-xj)
    Pcoi = 0;
    for i = 1:nloc
        for j = i+1:nloc
            Pcoi = Pcoi + D(i,j) * cos(x(i) - x(j));
        end
    end
    Pcoi = 2 * Pcoi; % because we only summed j>i

    % Build equality constraints: Pm - Pe + COI_terms = 0
    ceq = zeros(nloc,1);
    for i = 1:nloc
        % Electrical power Pe_i
        Pe_i = 0;
        for j = 1:nloc
            if i ~= j
                Pe_i = Pe_i + ( C(i,j) * sin(x(i) - x(j)) + D(i,j) * cos(x(i) - x(j)) );
            else
                Pe_i = Pe_i + E(i)^2 * real(Y1(i,i));
            end
        end
        % COI correction term:
        % (H_i / sum(H)) * ( Pcoi - sum(Pm - E.^2 .* real(diag(Y1))) )
        sys_term = sum(Pm - E.^2 .* real(diag(Y1))');
        coiterm = (H(i) / sum(H)) * ( Pcoi - sys_term );

        % ceq(i) = Pm(i) - Pe_i + COI_correction  -> solver sets ceq==0
        ceq(i) = Pm(i) - Pe_i + coiterm;
    end
end

%% ---------------- Utility ----------------
function xw = wrapToPi(x)
    % wrapToPi - wrap angle(s) to [-pi, pi]
    xw = mod(x + pi, 2*pi) - pi;
end
