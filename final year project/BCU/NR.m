function x = NR(f_func, J_func, x0, accuracy, max_itr)
% NR  Simple Newton-Raphson solver wrapper used by other scripts.
%
%   x = NR(f_func, J_func, x0)
%   x = NR(f_func, J_func, x0, accuracy, max_itr)
%
% f_func: handle returning f(x)
% J_func: handle returning J(x)
% x0     : initial guess vector
% accuracy : stopping tolerance (optional, default 1e-6)
% max_itr  : maximum iterations (optional, default 30)

if nargin < 4 || isempty(accuracy)
    accuracy = 1e-6;
end
if nargin < 5 || isempty(max_itr)
    max_itr = 30;
end

x = x0;
% preallocate history arrays for debugging (grow if needed)
nf = zeros(max_itr,1);
ndx = zeros(max_itr,1);

for iter = 1:max_itr
    J = J_func(x);          % Jacobian at current x
    f = f_func(x);          % residual at current x
    dx = - J \ f;           % Newton step

    nf(iter) = norm(f);
    ndx(iter) = norm(dx);

    % update solution
    x = x + dx;

    % stopping criteria (both residual and step magnitude)
    if (nf(iter) < accuracy) && (ndx(iter) < accuracy)
        % trim history arrays if desired
        break;
    end
end

if iter >= max_itr
    warning('NR:MaxIter', 'NR did not converge within max iterations');
    x = Inf;
    % optional debug plot
    try
        figure;
        subplot(2,1,1); plot(ndx(1:iter)); title('||dx||'); grid on;
        subplot(2,1,2); plot(nf(1:iter)); title('||f||'); grid on;
    catch
    end
end

end
