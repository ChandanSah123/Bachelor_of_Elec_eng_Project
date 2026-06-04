function Data1_Gen
%% 1. System Parameters (Extracted from your code)
Y = [1.1411-2.298*1i  0.1323+0.7035*1i  0.1854+1.0611*1i;
     0.1323+0.7035*1i 0.3810-2.0202*1i  0.1965+1.2031*1i;
     0.1854+1.0611*1i 0.1965+1.2031*1i  0.2723-2.3544*1i];

E = [1.0728 1.0775 1.0609];       % Internal Voltages (pu)
M = [23.64 6.4 3.01];             % Inertia Constants (2H)
Pm = [0.71 1.63 0.85];            % Mechanical Power (Assumed typical values for WSCC)
                                  % Note: Your code reads Pm from data, but usually Pm is constant.
                                  
ws = 1;                           % Synchronous speed (pu)
w_base = 377;                     % Base frequency (rad/s)

%% 2. Initial Conditions (Pre-fault / Disturbed State)
% We set initial angles slightly away from equilibrium to simulate a "swing"
delta0 = [0.1; 0.6; 0.4];         % Initial Angles (radians)
omega0 = [1.0; 1.0; 1.0];         % Initial Speeds (pu)

% Combine into state vector: [d1 w1 d2 w2 d3 w3]
x0 = [delta0(1); omega0(1); delta0(2); omega0(2); delta0(3); omega0(3)];

%% 3. Simulation (Solving the Swing Equation)
tspan = [0 2]; % Simulate for 2 seconds
options = odeset('RelTol',1e-5);
[t, x] = ode45(@(t,x) swing_equations(t, x, Pm, Y, E, M, ws, w_base), tspan, x0, options);

%% 4. Format Output into 'data1'
% Your code expects:
% Cols 1-3:   Angles (Degrees)
% Cols 4-6:   Mechanical Power (Pm)
% Cols 7-9:   (Unused/Electrical Power)
% Cols 10-12: Speed (pu)

num_steps = length(t);
data1 = zeros(num_steps, 12);

% --- Fill Columns ---
% 1. Angles: ODE outputs Radians, Code expects Degrees
data1(:, 1) = x(:, 1) * 180/pi; 
data1(:, 2) = x(:, 3) * 180/pi;
data1(:, 3) = x(:, 5) * 180/pi;

% 2. Mechanical Power (Constant)
data1(:, 4) = Pm(1);
data1(:, 5) = Pm(2);
data1(:, 6) = Pm(3);

% 3. Electrical Power (Optional - just filling space for indexing)
% We calculate this just so columns 7-9 aren't empty, though your code skips them.
for k = 1:num_steps
    delta = [x(k,1) x(k,3) x(k,5)];
    for i = 1:3
        Pe = 0;
        for j = 1:3
            Pe = Pe + E(i)*E(j)*abs(Y(i,j))*cos(angle(Y(i,j)) - (delta(i)-delta(j)));
        end
        data1(k, 6+i) = Pe; 
    end
end

% 4. Speed (pu)
% ODE output is the state variable.
data1(:, 10) = x(:, 2);
data1(:, 11) = x(:, 4);
data1(:, 12) = x(:, 6);

% Save to workspace
assignin('base', 'data1', data1);
disp('data1 generated successfully.');
disp(['Size: ', num2str(size(data1))]);

end

%% --- Helper Function: Swing Equation Dynamics ---
function dxdt = swing_equations(~, x, Pm, Y, E, M, ws, w_base)
    % Ensure x is a column vector
    x = x(:);

    dxdt = zeros(6,1);
    
    % Extract states
    d = x(1:2:end); % Angles: indices 1, 3, 5
    w = x(2:2:end); % Speeds: indices 2, 4, 6
    
    % Calculate Electrical Power (Pe)
    Pe = zeros(3,1);
    for i = 1:3
        sum_val = 0;
        for j = 1:3
            % Ybus Power Flow Equation
            % theta_ij = angle(Y(i,j))
            sum_val = sum_val + E(i)*E(j)*abs(Y(i,j))*cos(angle(Y(i,j)) - (d(i)-d(j)));
        end
        Pe(i) = sum_val;
    end
    
    % Damping (Assume zero or small)
    D = [0 0 0]; 
    
    % Differential Equations
    % 1. d(delta)/dt = (w - ws) * w_base
    dxdt(1:2:end) = (w - ws) * w_base;
    
    % 2. d(omega)/dt = (Pm - Pe - D(w-ws)) / M
    for k = 1:3
        idx = 2*k;
        dxdt(idx) = (Pm(k) - Pe(k) - D(k)*(w(k)-ws)) / M(k);
    end
end