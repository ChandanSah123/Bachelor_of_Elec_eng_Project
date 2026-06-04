%% BCU (Boundary Controlling u.e.p. Method)
tic; t_idx=1; BCU_idx=0;
max_BCU_sim=round(5/del_t_BCU);
theta_sim=zeros(num_gen,max_BCU_sim); theta_abs=zeros(1,max_BCU_sim);
theta_sim(:,1)=x_fault(idx_delta,PEBS_idx); theta_abs(1)=sum(abs(f_post_theta(theta_sim(:,1))));

while t_idx<max_BCU_sim && PEBS_idx~=0 && BCU_idx==0
    f=@(x) -x+theta_sim(:,t_idx)+del_t_BCU/2*(f_post_theta(theta_sim(:,t_idx))+f_post_theta(x)); % Trapezoidal rule
    J=@(x) -eye(num_gen)+del_t_BCU/2*J_post_theta(x);
    theta_sim(:,t_idx+1)=NR(f,J,theta_sim(:,t_idx));
    theta_abs(t_idx+1)=sum(abs(f_post_theta(theta_sim(:,t_idx+1))));
    t_idx=t_idx+1;
    if theta_abs(t_idx)>theta_abs(t_idx-1) && theta_abs(max(1,t_idx-2))>theta_abs(t_idx-1); BCU_idx=t_idx-2; end
end
t_bcu_sim=(0:t_idx-1)*del_t_BCU; theta_sim=theta_sim(:,1:t_idx); theta_abs=theta_abs(:,1:t_idx);

theta_u=theta_sim(:,BCU_idx);
x_eq_unstable=[theta_u; zeros(num_gen,1)];
Vcr_BCU=-P_post_kron'*(theta_u-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*theta_u)-cos(E_kron*x_eq_post(idx_delta)))+((abs(E_kron)*(theta_u-x_eq_post(idx_delta)))./(E_kron*(theta_u-x_eq_post(idx_delta))))'*diag(gij_post_kron)*(sin(E_kron*theta_u)-sin(E_kron*x_eq_post(idx_delta)))-V_pe0;
idx_tcr_BCU=find(V_total<Vcr_BCU); idx_tcr_BCU=idx_tcr_BCU(end);
tcr_BCU=del_t_fault*idx_tcr_BCU;
result_BCU=V_total(idx_tcl)<Vcr_BCU;
tcomp_BCU=toc; tcomp_BCU=tcomp_BCU+tcomp_PEBS;