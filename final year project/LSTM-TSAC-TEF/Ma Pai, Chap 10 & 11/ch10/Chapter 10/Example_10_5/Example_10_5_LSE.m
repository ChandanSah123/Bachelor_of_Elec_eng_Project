% Example 10.5   JHC June 5, 2017

% Linear State Estimator Code
% data is S (Bus 1) to F (Bus 2) to E (Bus 3)
% PMU data available for S and E

% line parameters (345 kV), per unit on 100 MVA base
% S to F
r_12 = 0.00003;
x_12 = 0.00042;
b_12 = 0.0;
% F to E
r_32 = 0.00251;
x_32 = 0.0346;
b_32 = 0.592/2;

% voltage base: 345 kV, current base = 100 MVA/345 kV/sqrt(3) = 167.3479
% raw PMU measurements at 45:34.0
% S voltage = 206368.1, phase = -15.912
% E (1) voltage = 204659.1, phase = -24.12
% E (2) voltage = 205259.6, phase = -23.936
% Line 10 (SF) current mag = 380.4, phase = 149.856
% Line FE1 current mag = 978.4, phase = 159.488, 
% PSE solution
% E voltage = 35668.89, phase = -101.0035
% S voltage = 778.0697, phase = -90
% F voltage = 778.0697, phase = -90
% Line 10 current mag = 380.3938, phase = 165.7656
% Line FE1 current mag = 869.27, phase 173.5104

% raw PMU data
% Bus 1 is S, Bus 2 is F, Bus 3 is E
jay = sqrt(-1);
V1p = 206368.1/345000*sqrt(3)*exp(jay*-(15.912)*pi/180);
%V1p = 206368.1/345000*sqrt(3)*exp(jay*-(14.65)*pi/180);  % need to add to
%current as well
V3p = 204659.1/345000*sqrt(3)*exp(jay*-24.12*pi/180);
I12p = 380.4/167.3479*exp(jay*149.856*pi/180);
% I32p = 978.4/167.3479*exp(jay*159.488*pi/180);
Is = 1
I32p = (978.4/Is)/167.3479*exp(jay*159.488*pi/180);



% A matrix set up
A = [ ...
    1 0 0 0 0 0 0 0 0 0; ...
    0 1 0 0 0 0 0 0 0 0; ...
    0 0 0 0 1 0 0 0 0 0; ...
    0 0 0 0 0 1 0 0 0 0; ...
    0 0 0 0 0 0 1 0 0 0; ...
    0 0 0 0 0 0 0 1 0 0; ...
    0 0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 0 1; ...
    (1-x_12*b_12) -r_12*b_12 -1 0 0 0 -r_12 x_12 0 0; ...
     r_12*b_12 (1-x_12*b_12) 0 -1 0 0 -x_12 -r_12 0 0; ...
     0 0 -1 0 (1-x_32*b_32) -r_32*b_32 0 0 -r_32 x_32; ... 
     0 0 0 -1 r_32*b_32 (1-x_32*b_32) 0 0 -x_32 -r_32  ...
    ];
b = [real(V1p) imag(V1p) real(V3p) imag(V3p) real(I12p) imag(I12p) ...
     real(I32p) imag(I32p) 0 0 0 0]';

% weighting matrix
W = eye(12,12);
W(5,5) = 1/2.27; W(6,6) = W(5,5); % W(9,9) = W(5,5); W(10,10) = W(5,5);
W(7,7) = 1/5.85; W(8,8) = W(7,7); % W(11,11) = W(7,7); W(12,12) = W(7,7);
W = eye(12,12);
W(9,9) = 100; W(10,10)=W(9,9); W(11,11)=W(9,9); W(12,12)=W(9,9);
% use high weights on ckt eqns
 
% x = (A'*W*A)\(A'*W*b);
x = (W*A)\(W*b);


% compare results
[b(1:8',1), x([1 2 5 6 7 8 9 10],1)];

% raw PMU data
V1 = 206368.1/345000*sqrt(3);
V1ph = -15.912;
V3 = 204659.1/345000*sqrt(3);
V3ph = -24.12;
I_12 = 380.4/167.3479;
I_12ph = 149.856;
I_32 = (978.4/Is)/167.3479;
I_32ph = 159.488;

% LSE results
V1_LSE = sqrt(x(1)^2+x(2)^2);
V1_ph_LSE = atan2(x(2),x(1))*180/pi;
V3_LSE = sqrt(x(5)^2+x(6)^2);
V3_ph_LSE = atan2(x(6),x(5))*180/pi;
I12_LSE = sqrt(x(7)^2+x(8)^2);
I12_ph_LSE = atan2(x(8),x(7))*180/pi;
I32_LSE = sqrt(x(9)^2+x(10)^2);
I32_ph_LSE = atan2(x(10),x(9))*180/pi;
V2_LSE = sqrt(x(3)^2+x(4)^2)
V2_ph_LSE = atan2(x(4),x(3))*180/pi

[V1 V1ph V1_LSE V1_ph_LSE]
[V3 V3ph V3_LSE V3_ph_LSE]
[I_12 I_12ph I12_LSE I12_ph_LSE]
[I_32 I_32ph I32_LSE I32_ph_LSE]

% I_12 calculated from LSE
I_12_PSE = (-V1_LSE*exp(jay*V1_ph_LSE*pi/180) + V2_LSE*exp(jay*V2_ph_LSE*pi/180))/(r_12+jay*x_12)


% results without current scaling:  Dec 2, 2015; use for Example 10.5

% Linear_PSE
% V2_LSE =    1.0373
% V2_ph_LSE =  -14.3246
% ans =    1.0361  -15.9120    1.0370  -14.3745
% ans =    1.0275  -24.1200    1.0312  -25.6393
% ans =    2.2731  149.8560    2.2731  149.8561
% ans =    5.8465  159.4880    5.8455  159.4868
% I_12_PSE =    1.9722 - 1.1430i












