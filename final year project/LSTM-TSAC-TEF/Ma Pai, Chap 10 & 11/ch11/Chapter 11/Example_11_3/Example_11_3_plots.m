% Example 11.3 JHC June 6, 2017

% make plots for Sauer&Pai, Example 11.3
% for both stable and unstable cases 
% Joe Chow Dec 16, 2015
% filename: Example_11_3_plots.m
load Example_11_3_stable
mac_ang_s = mac_ang; mac_spd_s = mac_spd; pelect_s = pelect;
load Example_11_3_unstable
mac_ang_u = mac_ang; mac_spd_u = mac_spd; pelect_u = pelect;

jay = sqrt(-1);
% plot power angle curve
% machine internal voltage
I = conj((bus_sol(1,4)+jay*bus_sol(1,5))/(bus_sol(1,2)*exp(jay*bus_sol(1,3)*pi/180)))
Ep = I*jay*0.245 + bus_sol(1,2)*exp(jay*bus_sol(1,3)*pi/180), abs(Ep), angle(Ep)*180/pi

theta = [0:1:180];
X_dp = 0.245; X_T = 0.15; X_e = 0.2;
P = abs(Ep)*1/(X_dp+X_T+X_e)*sin(theta*pi/180);
set_font_size
figure(1)
plot(theta,P), xlabel('Rotor angle (degrees)'), ylabel('Power transfer (pu)')
%title('P-delta curve for steady-state pre-fault condition')
hold on, plot([0 180],[0.9 0.9])
hold on, plot(mac_ang_s(1,:)*180/pi,pelect_s(1,:))
%title('P-delta curve for disturbance simulation')
%hold on, plot(atan(tan(mac_ang_u(1,:)))*180/pi,pelect_u(1,:)) %what a mess

figure(2)
plot(t,mac_ang_s(1,:)*180/pi), xlabel('Time (sec)'), ylabel('Rotor angle (degrees)')
%title('Machine angle response')
hold on, plot(t,mac_ang_u(1,:)*180/pi) 
axis([0 5 0 120])

figure(3)
plot(t,mac_spd_s(1,:)) 
xlabel('Time (sec)'), ylabel('Rotor speed (pu)')
%title('Machine speed response')
hold on, plot(t,mac_spd_u(1,:))
axis([0 5 0.98 1.03])
reset_font_size

%(ii)
% jay = sqrt(-1);
% I = conj((bus_sol(1,4)+jay*bus_sol(1,5))/(bus_sol(1,2)*exp(jay*bus_sol(1,3)*pi/180)))
% Ep = I*jay*0.245 + bus_sol(1,2)*exp(jay*bus_sol(1,3)*pi/180), abs(Ep), angle(Ep)*180/pi
% theta = [0:1:180];
% P = abs(Ep)*1/(X_dp+X_T+X_e)*sin(theta*pi/180);
% figure(1)
% plot(theta,P), xlabel('rotor angle in deg'), ylabel('power transfer in pu')
% hold on, plot([0 180],[0.9 0.9])

%(iii)
set_font_size
figure(4)
plot(mac_ang_s(1,:)*180/pi,mac_spd_s(1,:)-1)
xlabel('Rotor angle (degrees)'), ylabel('Rotor speed deviation (pu)')
%title('Phase plot of machine speed vs machine angle')
hold on, plot(mac_ang_u(1,:)*180/pi,mac_spd_u(1,:)-1)
axis([-10 110 -0.02 0.03])
reset_font_size

%(iv) post-fault equilibrium
X_epf = 0.4; %post-fault X_e
Xequi = X_dp+X_T+X_epf+0.01*991/1e9;
del_post = asin(0.9*Xequi/(abs(Ep)*1))*180/pi

% check VPE
del = [0:0.01:pi]';
VPE_del = -0.9*del - (abs(Ep)*1/Xequi)*cos(del) + (abs(Ep)*1/Xequi)*cosd(del_post);
%figure(10), plot(del,VPE_del)

%(v) plot of VPE, VKE, and VE

del_post = asin(0.9*Xequi/(abs(Ep)*1))
VPE_s = (-0.9*((mac_ang_s(1,:)-mac_ang_s(2,:))-del_post) - abs(Ep)/Xequi*(cos(mac_ang_s(1,:)-mac_ang_s(2,:))-cos(del_post))); 
VKE_s = 2.8756*2*pi*60*(mac_spd_s(1,:)-mac_spd_s(2,:)).*(mac_spd_s(1,:)-mac_spd_s(2,:));
VE_s  = VPE_s + VKE_s;
%BaseMVA = 991;
BaseMVA = 1;
set_font_size
figure(5), plot(t,BaseMVA*VPE_s,t,BaseMVA*VKE_s,t,BaseMVA*VE_s)
legend('VPE','VKE','VE')
%xlabel('Time (sec)'), ylabel('Energy functions (MW-sec)')
xlabel('Time (sec)'), ylabel('Energy functions ((pu-power)-sec)')
reset_font_size

