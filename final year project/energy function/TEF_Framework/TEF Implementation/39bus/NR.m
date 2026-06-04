function x = NR(f_func, J_func, x0, accuracy, max_itr)
% Robust Newton-Raphson Solver using Pseudo-Inverse
% Handles singular Jacobians common in IEEE 39-bus angle stability

if nargin < 4; accuracy = 1e-6; end
if nargin < 5; max_itr  = 50;   end

x = x0;
damping = 0.1; % Gentle damping prevents overshoot

for iter = 1:max_itr
    f = f_func(x);
    
    % Check convergence
    if norm(f, inf) < accuracy
        return;
    end
    
    J = J_func(x);
    
    % --- CRITICAL FIX FOR 39-BUS SYSTEM ---
    % The Jacobian is singular because angles are relative.
    % Standard division (\) fails. 'pinv' finds the minimum-norm solution.
    if rcond(J) < 1e-12
        dx = -pinv(J) * f;
    else
        dx = -J \ f;
    end
    % --------------------------------------
    
    x = x + damping * dx;
end

% Failure handling
disp(['NR did not converge. Final Residual: ' num2str(norm(f))]);
x = inf(size(x0)); 
end