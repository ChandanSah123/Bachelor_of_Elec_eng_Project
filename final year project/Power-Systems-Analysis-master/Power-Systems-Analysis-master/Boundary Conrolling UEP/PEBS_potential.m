%% PEBS (Potential Energy Boundary Surface)
tic; x_fault=x_eq_pre;
max_PEBS_sim=10000;
V_ke=zeros(1,max_PEBS_sim); V_p=zeros(1,max_PEBS_sim); V_d=zeros(1,max_PEBS_sim);
V_ke(1)=0.5*sum(M_gen.*x_fault(idx_omega,1).^2);
V_p(1)=-P_post_kron'*(x_fault(idx_delta,1)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_fault(idx_delta,1))-cos(E_kron*x_eq_post(idx_delta)));
V_d(1)=0;
V_pe0=-P_post_kron'*(x_eq_pre(idx_delta)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_eq_pre(idx_delta))-cos(E_kron*x_eq_post(idx_delta)))+0.5*(gij_post_kron.*(cos(E_kron*x_eq_pre(idx_delta))+cos(E_kron*x_eq_post(idx_delta))))'*(abs(E_kron)*x_eq_pre(idx_delta)-abs(E_kron)*x_eq_post(idx_delta));
V_total=V_ke+V_p+V_d;
PEBS=zeros(1,max_PEBS_sim);
PEBS(1)=f_post_theta(x_fault(idx_delta,1))'*(x_fault(idx_delta,1)-x_eq_post(idx_delta));
PEBS_idx=0; t_idx=1;
idx_tcl=round((t_cl+t_fault)/del_t_fault);
idx_tf=max(round(t_fault/del_t_fault),1);

while t_idx<max_PEBS_sim && (PEBS_idx==0 || t_idx<=idx_tcl)
    f=@(x) -x+x_fault(:,t_idx)+del_t_fault/2*(f_fault_kron(x_fault(:,t_idx))+f_fault_kron(x));
    J=@(x) -eye(num_var)+del_t_fault/2*J_fault_kron(x);
    x_fault(:,t_idx+1)=NR(f,J,x_fault(:,t_idx));
    V_ke(t_idx+1)=0.5*sum(M_gen.*x_fault(idx_omega,t_idx+1).^2);
    V_p(t_idx+1)=-P_post_kron'*(x_fault(idx_delta,t_idx+1)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_fault(idx_delta,t_idx+1))-cos(E_kron*x_eq_post(idx_delta)));
    V_d(t_idx+1)=V_d(t_idx)+0.5*(gij_post_kron.*(cos(E_kron*x_fault(idx_delta,t_idx+1))+cos(E_kron*x_fault(idx_delta,t_idx))))'*(abs(E_kron)*x_fault(idx_delta,t_idx+1)-abs(E_kron)*x_fault(idx_delta,t_idx));
    V_total(t_idx+1)=V_ke(t_idx+1)+V_p(t_idx+1)+V_d(t_idx+1);
    PEBS(t_idx+1)=f_post_theta(x_fault(idx_delta,t_idx+1))'*(x_fault(idx_delta,t_idx+1)-x_eq_post(idx_delta));
    if PEBS(t_idx)<0 && PEBS(t_idx+1)>=0 && PEBS_idx==0; PEBS_idx=t_idx; end
    t_idx=t_idx+1;
end
t_fault_sim=(0:t_idx-1)*del_t_fault; V_ke=V_ke(:,1:t_idx); V_p=V_p(:,1:t_idx); V_d=V_d(:,1:t_idx); V_total=V_total(:,1:t_idx); PEBS=PEBS(:,1:t_idx);

Vcr_PEBS=V_p(PEBS_idx)+V_d(PEBS_idx)-V_pe0; 
idx_tcr_PEBS=find(V_total<Vcr_PEBS); idx_tcr_PEBS=idx_tcr_PEBS(end);
tcr_PEBS=del_t_fault*idx_tcr_PEBS;
result_PEBS=V_total(idx_tcl)<Vcr_PEBS;
tcomp_PEBS=toc;