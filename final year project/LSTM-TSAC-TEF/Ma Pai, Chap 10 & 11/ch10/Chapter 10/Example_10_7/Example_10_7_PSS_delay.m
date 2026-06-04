% m file for Example 10.7 - PSS design with time lag

% first check the matrices in Example 8.7
% no PSS
A1 = [ -0.3511 -0.236   0  0.104; ...
        0       0     377  0 ; ...
       -0.1678 -0.144   0  0 ; ...
       -714.4 -10       0  -5];
eig(A1)

%   -2.5881 + 8.4950i
%   -2.5881 - 8.4950i
%   -0.0875 + 7.1110i
%   -0.0875 - 7.1110i

A2 = [ -0.3511 -0.236   0  0.104  0; ...
        0       0     377  0      0; ...
       -0.1678 -0.144   0  0      0; ...
       -714.4  -10      0  -5    2000; ...
       -0.42   -0.36    5  0     -10];
eig(A2)

%   -1.6314 + 8.5504i
%   -1.6314 - 8.5504i
%   -0.8612 + 7.0742i
%   -0.8612 - 7.0742i
%  -10.3660          

% set up system to repeat root locus design of Example 8.7
B1 = [0; 0; 0; 2000];
C1 = [0 0 1 0];
Gep = ss(A1,B1,C1,0);
Gs = tf([0.5 1],[0.1 1]); 
set_font_size
figure(1), rlocus(-Gep*Gs)
axis([-3 1 -2 12])
reset_font_size

% add time delay
Td = 0.1
pade_tf = tf([-Td/2 1],[Td/2 1])
set_font_size
figure(2), rlocus(-pade_tf*Gep*Gs)
axis([-3 1 -2 12])
reset_font_size

% phase lag at 7 rad/s
ph_lag = 7*.1*180/pi   % 40 degrees

% 40 degrees phase lead => pole-zero ration of 5
% 1/sqrt(0.5*0.1)
T3 = 1/(7/sqrt(5)); T4 = 1/(7*sqrt(5));  %  T3 = 0.3194, T4 = 0.0639
Gs2 = tf([T3 1],[T4 1]); 
set_font_size
figure(3), rlocus(-pade_tf*Gep*Gs*Gs2)
axis([-3 1 -2 12])
reset_font_size

CL_syst = feedback(Gep*pade_tf,-Gs*Gs2*0.2)
pole(CL_syst)

%  -19.0066 + 4.4861i
%  -19.0066 - 4.4861i
%   -1.4388 + 8.6758i
%   -1.4388 - 8.6758i
%   -0.7702 + 6.9848i
%   -0.7702 - 6.9848i
%   -8.5723          







