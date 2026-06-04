%% Set up variables
run(['dyn' int2str(sys_case)]) %1. Get info from data file/load system data
slack_bus=SW.con(1);           %2.Extract system sizes
num_bus=size(Bus.con,1);
num_line=size(Line.con,1);
num_gen=size(Syn.con,1);
num_load=num_bus-num_gen;
%3.define indices
eye_bus=eye(num_bus);
idx_gen=Syn.con(:,1);
idx_load=setdiff(1:num_bus,idx_gen);
idx_delta=1:num_gen;
idx_omega=num_gen+1:2*num_gen;
num_var=2*num_gen;
%4.Generator parameters
M_gen=Syn.con(:,18)./(2*pi*Syn.con(:,4));
M_T=sum(M_gen);
Syn.con(:,19)=2*ones(size(Syn.con(:,19))); % Add damping
D_gen=Syn.con(:,19)./(2*pi*Syn.con(:,4));  %damping normalizing
v_gen=[SW.con(:,4); PV.con(:,5)];          %generator voltage setpoing
xd_p=Syn.con(:,9);                          %transient reactance Xd'

%5.Build the raw network admittance matrix Y
 %This constructs the network Ybus using primitive incidence-based stamping:
line_frto=Line.con(:,1:2);
fault_bus=line_frto(fault_line,fault_frto_bus);
Z_line=Line.con(:,8)+1i*Line.con(:,9);
E=eye_bus(line_frto(:,1),:)-eye_bus(line_frto(:,2),:);
Y=E'*diag(Z_line.^-1)*E; 

%6.Generator injections and load powers
Pgen=zeros(num_bus,1);
Pgen(PV.con(:,1))=PV.con(:,4); %sets generator real powers (PV buses)
Sload=zeros(num_bus,1); 
Sload(PQ.con(:,1))=(PQ.con(:,4)+1i*PQ.con(:,5)); %sets load complex powesrs(PQ)

%7.Convert to constant impedence load/ structure preserving load
x_eq=NR_ss(Y,Pgen-Sload,idx_load,v_gen,slack_bus); %solves load flow

V_eq=x_eq(num_bus+1:end).*(cos(x_eq(1:num_bus))+1i*sin(x_eq(1:num_bus)));
I_eq=Y*V_eq;
S_inj=V_eq.*conj(I_eq); %compute complex bus voltages and currents.

y_load=conj(Sload)./V_eq.^2; %converts PQ load into constand admittance Yload


%8.Build pre_fault and Post_fault network Y
YN_pre=E'*diag(Z_line.^-1)*E+diag(y_load);  %add load as admittances
YN_fault=YN_pre;
Zf_line=Z_line;
Zf_line(fault_line)=inf;  % infinite impedance for fault line
YN_post=E'*diag(Zf_line.^-1)*E+diag(y_load); %one line removed

%9.Include stator impedence for network reduction/include generator
%transient reactance Xd'
Y_pre=zeros(num_bus+num_gen);
Y_pre([1:num_gen,num_gen+idx_gen'],[1:num_gen,num_gen+idx_gen'])=[diag((1i*xd_p).^-1) diag(-(1i*xd_p).^-1); diag(-(1i*xd_p).^-1) diag((1i*xd_p).^-1)];
%This adds
%| 1/jXd'   -1/jXd' |
%| -1/jXd'   1/jXd' |

Y_fault=Y_pre; 
Y_post=Y_pre;
Y_pre(num_gen+1:end,num_gen+1:end)=Y_pre(num_gen+1:end,num_gen+1:end)+YN_pre;
Y_fault(num_gen+1:end,num_gen+1:end)=Y_fault(num_gen+1:end,num_gen+1:end)+YN_fault;
Y_fault(num_gen+fault_bus,:)=[]; Y_fault(:,num_gen+fault_bus)=[];
Y_post(num_gen+1:end,num_gen+1:end)=Y_post(num_gen+1:end,num_gen+1:end)+YN_post;

%10. Pre-contingency Equilibrium/ calculate pre and post-fault equilibrium
x_eq_pre=NR_ss(YN_pre,Pgen,idx_load,v_gen,slack_bus);
V_eq_pre=x_eq_pre(num_bus+1:end).*(cos(x_eq_pre(1:num_bus))+1i*sin(x_eq_pre(1:num_bus)));
I_eq_pre=YN_pre*V_eq_pre;
Pgen_pre=real(V_eq_pre(idx_gen).*conj(I_eq_pre(idx_gen)));

%internal EMF computation
Eeq_pre=abs(V_eq_pre(idx_gen)+1i*xd_p.*I_eq_pre(idx_gen));
delta_eq_pre=angle(V_eq_pre(idx_gen)+1i*xd_p.*I_eq_pre(idx_gen));

%constructing COI_reference states
x_eq_pre=[delta_eq_pre-M_gen'*delta_eq_pre/M_T; zeros(num_gen,1)];

% Post-contingency Equilibrium
x_eq_post=NR_ss(YN_post,Pgen,idx_load,v_gen,slack_bus);
V_eq_post=x_eq_post(num_bus+1:end).*(cos(x_eq_post(1:num_bus))+1i*sin(x_eq_post(1:num_bus)));
I_eq_post=YN_post*V_eq_post;
Pgen_post=real(V_eq_post(idx_gen).*conj(I_eq_post(idx_gen)));
Eeq_post=abs(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
delta_eq_post=angle(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
x_eq_post=[delta_eq_post-M_gen'*delta_eq_post/M_T; zeros(num_gen,1)];

% 11.Apply Kron Reduction
Y_pre_kron=Y_pre(1:num_gen,1:num_gen)-Y_pre(1:num_gen,num_gen+1:end)*(Y_pre(num_gen+1:end,num_gen+1:end)\Y_pre(num_gen+1:end,1:num_gen));
Y_fault_kron=Y_fault(1:num_gen,1:num_gen)-Y_fault(1:num_gen,num_gen+1:end)*(Y_fault(num_gen+1:end,num_gen+1:end)\Y_fault(num_gen+1:end,1:num_gen));
Y_post_kron=Y_post(1:num_gen,1:num_gen)-Y_post(1:num_gen,num_gen+1:end)*(Y_post(num_gen+1:end,num_gen+1:end)\Y_post(num_gen+1:end,1:num_gen));

%12. Build edge_based E_kron
edge_kron=nchoosek(1:num_gen,2);
E_kron=zeros(size(edge_kron,1),num_gen);

%13.Build Bij and Gij
for i=1:size(E_kron,1)
    E_kron(i,edge_kron(i,1))=1;  E_kron(i,edge_kron(i,2))=-1;
    gij_pre_kron(i,1)=real(Y_pre_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_pre(E_kron(i,:)==1)*Eeq_pre(E_kron(i,:)==-1);
    gij_fault_kron(i,1)=real(Y_fault_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
    gij_post_kron(i,1)=real(Y_post_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
    bij_pre_kron(i,1)=imag(Y_pre_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_pre(E_kron(i,:)==1)*Eeq_pre(E_kron(i,:)==-1);
    bij_fault_kron(i,1)=imag(Y_fault_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
    bij_post_kron(i,1)=imag(Y_post_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
end
%14Build P_i terms
P_pre_kron=Pgen_pre-real(diag(Y_pre_kron)).*Eeq_pre.^2;
P_fault_kron=Pgen_post-real(diag(Y_fault_kron)).*Eeq_post.^2;
P_post_kron=Pgen_post-real(diag(Y_post_kron)).*Eeq_post.^2;