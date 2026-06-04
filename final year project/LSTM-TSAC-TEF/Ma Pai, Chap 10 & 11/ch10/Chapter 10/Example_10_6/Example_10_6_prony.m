% Example 10.6 

% Aug 10 is oscillation
[num10,txt10,raw10] = xlsread('August10_Joe_59-74.csv','A2:D302');
% Channel 1 (column) is time, channel 2 is Malin voltage, channel 3 is COI
% (Cal-Oregon interface), channel 4 is frequency
t = num10(:,1); V_Malin = num10(:,2); P_COI = num10(:,3); f = num10(:,4);
set_font_size
figure, plot(t,P_COI)
xlabel('Time (sec)'), ylabel('Active Power Flow (MW)')
axis([0 90 3400 5000])
reset_font_size

% Remove average flow
P_MW2 = P_COI(:,1) - 4200; % ends at 74 seconds


% the following are the results used in the chapter
n = 100  % 50 and 30 are good; 5 is no good
ts = 0.05;
%figure, plot(P_MW1)
[B,A] = prony(P_MW2, n-1, n);   % MATLAB prony function
syst = tf(B,A,ts)
[pole(syst) abs(pole(syst))]
td2 = 300*0.05;
[y,td] = impulse(syst,td2);
set_font_size
figure, plot(td+59,P_MW2,'k-',td+59,y*ts,'k--')
xlabel('Time (sec)'), ylabel('Active Power Flow (MW)')
legend('Raw data','Prony model')
axis([59 74 -200 250])
reset_font_size


z = 1.0004 + j*0.0714  % pole from n = 30, unstable
s = log(z)/ts
% z = 1.000391058749736e+00 + 7.136833526857518e-02i
% |z| = 1.00293355
% s = 0.0588 + j*1.4250  ==> 0.2268 Hz, damping ratio = -0.0412 (unstable)
syst_ss = canon(syst,'modal'); % mode is 52,53
syst_small = ss(syst_ss.a(52:53,52:53),syst_ss.b(52:53,1), ...
                syst_ss.c(1,52:53),0,ts)
y2 = impulse(syst_small,td2);
set_font_size
figure, plot(td+59,P_MW2,'k-',td+59,y2*ts,'k--')
xlabel('Time (sec)'), ylabel('Active Power Flow (MW)')
legend('Raw data','Second-order model')
axis([59 74 -200 250])
reset_font_size

% abs = 1
% if abc == 0
% syst_small =
%   a = 
%              x1        x2
%    x1         1   0.07137
%    x2  -0.07137         1
%  
%   b = 
%           u1
%    x1  21.85
%    x2  -25.3
%  
%   c = 
%            x1      x2
%    y1  -1.581  -1.861
%  
%   d = 
%        u1
%    y1   0
%  
% Sample time: 0.05 seconds
% Discrete-time state-space model.
% 
% tf(syst_small)
% 
% ans =
%  
%      12.54 z - 6.786   (12.5367  -6.7863)
%   ---------------------
%   z^2 - 2.001 z + 1.006  (-2.0008  1.0059)
%  
% Sample time: 0.05 seconds
% Discrete-time transfer function.
% end
