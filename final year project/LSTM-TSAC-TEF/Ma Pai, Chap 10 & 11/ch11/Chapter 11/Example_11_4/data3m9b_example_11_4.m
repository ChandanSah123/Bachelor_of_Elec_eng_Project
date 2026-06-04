% A 3-machine 9-bus system from Chow's slow coherency book p. 70
% The bus numbers are based on p. 165 of Sauer-Pai book
% data3m9b_example_11_4.m
% system base is 100 MVA

% bus data format
% bus: number, voltage(pu), angle(degree), p_gen(pu), q_gen(pu),
%      p_load(pu), q_load(pu),G shunt,B shunt, bus_type
%      bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)

bus = [ ...
%   #  V      phase  P_gen  Q_gen P_ld  Q_ld  G_sh  B_sh type
    1 1.04    0.00   0.00   0.00  0.00  0.00  0.00  0.00   1;
	2 1.025   0.00   1.63   0.00  0.00  0.00  0.00  0.00   2;
	3 1.025   0.00   0.85   0.00  0.00  0.00  0.00  0.00   2;
    4 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00   3;
    5 1.00    0.00   0.00   0.00  1.25  0.50  0.00  0.00   3;
	7 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00   3;
	8 1.00    0.00   0.00   0.00  1.00  0.35  0.00  0.00   3;
	9 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00   3];

% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio

line = [ ...
%   from to  R      X      B      tap phase
    1    4   0.0    0.0576 0.     0. 0. ;
	4    6   0.017  0.092  0.158  0. 0. ;
	6    9   0.039  0.17   0.358  0. 0. ;
	3    9   0.0    0.0586 0.     0. 0. ;
	9    8   0.0119 0.1008 0.209  0. 0. ;
	8    7   0.0085 0.072  0.149  0. 0. ;
	7    2   0.0    0.0625 0.     0. 0. ;
	7    5   0.032  0.161  0.306  0. 0. ;
	5    4   0.01   0.085  0.176  0. 0. ];


% Machine data format
% Machine data format
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant
%                T"_do(sec),
%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant
%                T"_qo(sec),
%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. bus number
%
% note: all the following machines use electro-mechanical model
mac_con = [ ...
% G# B# MVA xl ra xd      xd'     xd" Tdo' Tdo" xq     xq'    xq" Tq'   Tq" H     D D  B#           
% data from Pai book
%   1  1  100 0  0  0.146   0.0608  0   8.96 0    0.0969 0.0969 0   0.31  0   23.64 0 0  1;
%   2  2  100 0  0  0.8958  0.1198  0   6.0  0    0.8645 0.1969 0   0.535 0   6.4   0 0  2;
%   3  3  100 0  0  1.3125  0.1813  0   5.89 0    1.2578 0.25   0   0.6   0   3.01  0 0  3];
% change xq' to be the same as xd' (nelgecting transient saliency
  1  1  100 0  0  0.146   0.0608  0   8.96 0    0.0969 0.0608 0   0.31  0   23.64 0 0  1;
  2  2  100 0  0  0.8958  0.1198  0   6.0  0    0.8645 0.1198 0   0.535 0   6.4   0 0  2;
  3  3  100 0  0  1.3125  0.1813  0   5.89 0    1.2578 0.1813 0   0.6   0   3.01  0 0  3];


%
%       Exciter data
%
%       format
%column data
% 1     exciter type(1 for DC1,2 for DC2)
% 2     machine number
% 3     input filter time constant
% 4     voltage regulator gain K_A
% 5     voltage regulator time constant T_A(sec)
% 6     voltage regulator time constant T_B(sec)
% 7     voltage regulator time constant T_C(sec)
% 8     maximum voltage regulator output VR_max 
% 9     minimum voltage regulator output VR_min
% 10    exciter constat K_E
% 11    exciter time constant T_E
% 12    E_1
% 13    saturation function S_E(E_1)
% 14    E_2
% 15    saturation function S_E(E_2)
% 16    stabilizer gain K_F
% 17    stabilizer time constant(T_F)

exc_con=[...
% 1 is DC1A
% type # T_R K_A  T_A  T_B T_C  VRma VRmi  K_E    T_E   E_1 S_E1   E_2 S_E2   K_F    T_F 
  1    1  0 20.0  0.2  0   0   7.0   -7.0  1.0    0.314 0.0 0.0039 0.0 1.555 0.063  0.35;
  1    2  0 20.0  0.2  0   0   7.0   -7.0  1.0    0.314 0.0 0.0039 0.0 1.555 0.063  0.35;
  1    3  0 20.0  0.2  0   0   7.0   -7.0  1.0    0.314 0.0 0.0039 0.0 1.555 0.063  0.35];


% loads - constant power loads
load_con = [ ...
%             CP CP CI CI 
            5 1 1 0 0; ...
            6 1 1 0 0; ...
            8 1 1 0 0];

%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  - 0 three phase
%                            - 1 line to ground
%                            - 2 line-to-line to ground
%                            - 3 line-to-line
%                            - 4 loss of line with no fault
%                            - 5 loss of load at bus
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
%
%
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)

sw_con = [...
0    0    0    0    0    0    0.005;  % sets intitial time step
0.1  4    5    0    0    0    0.0025; % apply three phase fault at bus 4, on line 4-5
0.15 0    0    0    0    0    0.0025; % clear fault at bus 4
0.20 0    0    0    0    0    0.0025; % clear remote end
0.50 0    0    0    0    0    0.005;  % increase time step 
1.0  0    0    0    0    0    0.01;   % increase time step
5.0  0    0    0    0    0    0];     % end simulation


% some useful command

% run power flow 
% [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,30,1,'y',1);

% To see A matrix with small entries removed 
% a_mat_clean = small20(a_mat,1e-6);

% Build A matrix from c_curd, c_curq - compare to Pai's results 

% xx = [0.146   0.0608      0.0969 0.0969; ...
%       0.8958  0.1198      0.8645 0.1969; ... 
%       1.3125  0.1813      1.2578 0.25 ]  
% 
% H = [23.64  6.4  3.01]
% id = [0.3026  1.2901  0.5615]
% iq = [0.6712  0.9320  0.6194]
% edp = [0.0242  0.6941  0.6668]
% eqp = [1.0564  0.7782  0.7679]
% Tdop = [8.96  6.0    5.89]
% Tqop = [0.31  0.535  0.6]
% 
% BB1 = [0 0 0 0 0 0; ...
%        (iq(1)*(xx(1,2)-xx(1,4))-edp(1))/(2*H(1)) 0 0 ...
%        (id(1)*(xx(1,2)-xx(1,4))-eqp(1))/(2*H(1)) 0 0; ...
%        -(xx(1,1)-xx(1,2))/Tdop(1) 0 0 0 0 0; ...
%        0 0 0 (xx(1,3)-xx(1,4))/Tqop(1) 0 0; ...
%        0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; 
% 
% BB2 = [0 0 0 0 0 0; ...
%        0 (iq(2)*(xx(2,2)-xx(2,4))-edp(2))/(2*H(2)) 0 ...
%        0 (id(2)*(xx(2,2)-xx(2,4))-eqp(2))/(2*H(2)) 0 ; ...
%        0 -(xx(2,1)-xx(2,2))/Tdop(2) 0 0 0 0 ; ...
%        0 0 0 0 (xx(2,3)-xx(2,4))/Tqop(2) 0 ; ...
%        0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; 
%    
% BB3 = [0 0 0 0 0 0; ...
%        0 0 (iq(3)*(xx(3,2)-xx(3,4))-edp(3))/(2*H(3)) ...
%        0 0 (id(3)*(xx(3,2)-xx(3,4))-eqp(3))/(2*H(3)) ; ...
%        0 0 -(xx(3,1)-xx(3,2))/Tdop(3) 0 0 0 ; ...
%        0 0 0 0 0 (xx(3,3)-xx(3,4))/Tqop(3) ; ...
%        0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; 
%    
% BB = [BB1; BB2; BB3]; 
% 
% BVT = [zeros(4,3); ...
%        -100 0 0; ...
%        0 0 0; 0 0 0; ...
%        zeros(4,3); ...
%        0 -100 0; ...
%        0 0 0; 0 0 0; ...
%        zeros(4,3); ...
%        0 0 -100; ...
%        0 0 0; 0 0 0 ]; 
% 
% BB*[c_curd; c_curq] + BVT*c_v(1:3,:)     
    

