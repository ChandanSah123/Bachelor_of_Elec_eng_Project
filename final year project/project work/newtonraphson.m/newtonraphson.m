function [sol, error] = newtonraphson(F, x0, maxiter)
error1 = 1e8;
x = x0;
iter = 0;
if nargin < 3
    maxiter = 1e3;
end
error = zeros(maxiter, 1);

sol = F(x0); % compute test solution
nvar = length(sol); % compute size of the system
h = 1e-4 .* ones(nvar, 1); % size of step for derivative computation

while error1 > 1e-12
    iter = iter+1; % update iteration
    f = F(x); % evaluate function at current point
    
    J = jacobiannum(F, x, h); % compute the Jacobian
    y = -J\f;  % solve the linear equations
    x = x + y; % move the solution
    
    % calculate errors
    error1 = sqrt(sum(y.^2));
    error(iter) = sqrt(sum(f.^2));
    
    % break computation if maximum number of iterations is exceeded
    if iter == maxiter
        warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
        break;
    end 
end

% return solution and error for each iteration
sol = x;
error = error(1:iter);

end

function J = jacobiannum(F,x,h)
% Computes the Jacobian matrix of the function F at the point x, where h is
% the step size to take on each dimension (has to be small enough for
% precision). Note that the result of f and vectors x and h must be column
% vectors of the same dimensions.

J = (F(repmat(x,size(x'))+diag(h))-F(repmat(x,size(x'))))./h';

end