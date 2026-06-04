function x = NR_ss(Y, S_inj, idx_pq, v_PV, slack_bus)
% NR_ss  Newton-Raphson power flow solver (steady-state)
%
% x = NR_ss(Y, S_inj, idx_pq, v_PV, slack_bus)
%
% Inputs:
%   Y       - system admittance matrix (n x n)
%   S_inj   - complex injection vector (Sgen - Sload) (n x 1)
%   idx_pq  - indices of PQ buses (vector)
%   v_PV    - vector of voltage setpoints for PV and slack buses
%   slack_bus - index of slack bus (scalar)
%
% Output:
%   x = [theta; Vmag] where theta (rad) is bus angle vector and Vmag is voltage magnitudes

max_iter = 25;
num_bus = size(Y,1);

% Initialize voltage magnitudes (close to 1 p.u.) and angles
v_mag = 1 + 0.01 * rand(num_bus,1);         % initial magnitudes
v_mag(setdiff(1:num_bus, idx_pq)) = v_PV;  % set PV/slack voltage magnitudes

theta = 0.01 * rand(num_bus,1);             % initial angles
theta(slack_bus) = 0;                       % fix slack angle to zero

% set up index sets (exclude slack from angle unknowns)
idx_nslack = [1:slack_bus-1, slack_bus+1:num_bus];

for iter = 1:max_iter
    % complex bus voltages
    v_cpx = v_mag .* (cos(theta) + 1i * sin(theta));
    % complex power imbalance S = V * conj(Y*V) - S_inj
    S_bal = v_cpx .* conj(Y * v_cpx) - S_inj;

    % Form mismatch vector f = [P_mismatch (exclude slack); Q_mismatch (PQ buses only)]
    f = [ real(S_bal(idx_nslack)); imag(S_bal(idx_pq)) ];

    % Build Jacobian blocks following standard power-flow derivation
    % Note: these blocks are computed in compact vectorized form
    J1 = real( diag(v_cpx) * conj( Y .* (ones(num_bus,1) * (1i * v_cpx).') ) + diag(1i * v_cpx) * diag( conj(Y * v_cpx) ) );
    J2 = real( diag(v_cpx) * conj( Y .* (ones(num_bus,1) * (v_cpx ./ v_mag).') ) + diag(v_cpx ./ v_mag) * diag( conj(Y * v_cpx) ) );
    J3 = imag( diag(v_cpx) * conj( Y .* (ones(num_bus,1) * (1i * v_cpx).') ) + diag(1i * v_cpx) * diag( conj(Y * v_cpx) ) );
    J4 = imag( diag(v_cpx) * conj( Y .* (ones(num_bus,1) * (v_cpx ./ v_mag).') ) + diag(v_cpx ./ v_mag) * diag( conj(Y * v_cpx) ) );

    % Assemble reduced Jacobian for unknowns [theta(no-slack); Vmag(PQ)]
    J = [ J1(idx_nslack, idx_nslack), J2(idx_nslack, idx_pq);
          J3(idx_pq, idx_nslack),     J4(idx_pq, idx_pq) ];

    % Solve linear system
    dx = - J \ f;

    % Update angles and voltage magnitudes
    theta(idx_nslack) = theta(idx_nslack) + dx(1:length(idx_nslack));
    v_mag(idx_pq)      = v_mag(idx_pq)      + dx(length(idx_nslack)+1:end);

    % convergence check
    if (norm(f) < 1e-8) && (norm(dx) < 1e-8)
        break;
    end
end

% pack solution as [theta; Vmag]
x = [theta; v_mag];

if iter == max_iter
    warning('NR_ss:MaxIter', 'Newton Raphson did not fully converge');
    x = Inf;
end

end
