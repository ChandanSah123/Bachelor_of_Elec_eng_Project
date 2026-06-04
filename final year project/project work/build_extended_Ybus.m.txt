%%---------Extented Ybus for transient Stability ----------%%
function Yext = build_extended_Ybus(Ynet, gen_bus, Xd_val)
% Ynet    : N x N base network Ybus (lines/transformers/shunts)
% gen_bus : ng x 1 terminal bus indices (1..N)
% Xd_val  : ng x 1 vector of reactances in p.u. (choose Xd' for transient,
%           or Xd'' for fault)
%
% returns Yext of size (N+ng)x(N+ng) with internal nodes appended (N+1..N+ng)
N  = size(Ynet,1);
ng = length(gen_bus);
Yext = complex(zeros(N+ng, N+ng));
Yext(1:N,1:N) = Ynet;

for k = 1:ng
    tb = gen_bus(k);
    in = N + k;
    Yg = 1/(1j * Xd_val(k));    % admittance of generator reactance
    % off-diagonal
    Yext(tb, in) = Yext(tb, in) - Yg;
    Yext(in, tb) = Yext(tb, in);
    % diagonals
    Yext(tb, tb) = Yext(tb, tb) + Yg;
    Yext(in, in) = Yext(in, in) + Yg;
end
end
