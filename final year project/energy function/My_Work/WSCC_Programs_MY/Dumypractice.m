% Example (3-machine)
H = [23.64 6.4 3.01];
E = [1.0728 1.0775 1.0609];
%Pm = [0.7195 1.63 0.85];
Pm=[0.275542318820953	0.525806248188019	0.303571462631226];

% Example Yint_post (you must supply actual complex matrix)
% For demonstration I'll use a placeholder symmetric matrix (replace with your Y1)
Y1 = Yint_post;

% initial guess
x0=data(100,1:3);

[x_coi, x_raw, exitflag, lambda, resid] = solveSEP(Y1, Pm, E, H, x0)

% x_coi is the SEP angles in COI frame
