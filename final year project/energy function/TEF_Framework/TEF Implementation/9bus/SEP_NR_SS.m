clc;clear;close all;
load('Y_all.mat')
load("data1.mat")
num_bus=9;
num_gen=3;
idx_gen=[1 2 3];
idx_load=setdiff(1:num_bus,idx_gen);
v_gen=[1.04;1.025;1.025];
slack_bus=1;
H= [23.64;6.4;3.01];
M_gen=H.*(1/(pi*60));
xd_p=[0.0608;0.1198;0.1813];
M_T=sum(M_gen);
%Pm=data(10,4:6)';
Pm = [0;1.63;0.85;0;0;0;0;0;0];
Pgen=Pm;
YN_post=Y_post;
x_eq_post=NR_ss(YN_post,Pgen,idx_load,v_gen,slack_bus);
V_eq_post=x_eq_post(num_bus+1:end).*(cos(x_eq_post(1:num_bus))+1i*sin(x_eq_post(1:num_bus)));
I_eq_post=YN_post*V_eq_post;
Pgen_post=real(V_eq_post(idx_gen).*conj(I_eq_post(idx_gen)));
Eeq_post=abs(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
delta_eq_post=angle(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
x_eq_post=[delta_eq_post-M_gen'*delta_eq_post/M_T; zeros(num_gen,1)];

