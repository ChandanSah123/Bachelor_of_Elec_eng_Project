function dy = integrateth(t, y)
% -------------------------------------------------------------
% Three-machine COI-reduced swing equations.
% y(i) = rotor angle δ_i (in radians)
% dy(i) = dδ_i/dt = ω_i - ω_COI (implicitly embedded)
% -------------------------------------------------------------

% Generator internal voltages
E = [1.0728 1.0775 1.0609];

% Inertia constants
H = [23.64 6.4 3.01];

% Mechanical power input
Pm = [0.7195 1.63 0.85];

% Admittance matrix (internal bus) — no outage case
Y1 = [
    0.8450 - 2.9880i, 0.2870 + 1.5130i, 0.2100 + 1.2260i;
    0.2870 + 1.5130i, 0.4200 - 2.7240i, 0.2130 + 1.0880i;
    0.2100 + 1.2260i, 0.2130 + 1.0880i, 0.2770 - 2.3680i
];

g = 3;

% -------------------------------------------------------------
% Precompute coupling coefficients
% C(i,j) = E_iE_j imag(Y_ij)
% D(i,j) = E_iE_j real(Y_ij)
% -------------------------------------------------------------
C = zeros(g);
D = zeros(g);

for i = 1:g
    for j = 1:g
        C(i,j) = E(i)*E(j)*imag(Y1(i,j));
        D(i,j) = E(i)*E(j)*real(Y1(i,j));
    end
end

% -------------------------------------------------------------
% COI term calculations
% -------------------------------------------------------------
Pe_total = 0;
for i = 1:g
    Pe_total = Pe_total + (Pm(i) - real(Y1(i,i))*E(i)^2);
end

% Double-sum of D(i,j)cos(δ_i - δ_j)
Dcos = D(1,2)*cos(y(1) - y(2)) + ...
       D(1,3)*cos(y(1) - y(3)) + ...
       D(2,3)*cos(y(2) - y(3));

COI_term = Pe_total - 2*Dcos;

% -------------------------------------------------------------
% Swing equations in COI frame
% dy(i) = dδ_i/dt = f(δ)
% -------------------------------------------------------------
dy = zeros(3,1);

for i = 1:g
    % Electrical output Pe_i
    Pe_i = real(Y1(i,i))*E(i)^2;
    for j = 1:g
        if i ~= j
            Pe_i = Pe_i + C(i,j)*sin(y(i)-y(j)) + D(i,j)*cos(y(i)-y(j));
        end
    end

    % COI correction term
    COI_i = (H(i)/sum(H)) * COI_term;

    dy(i) = Pm(i) - Pe_i - COI_i;
end
end
