%% Direct Simulation
tic; num_sim_idx=ceil((t_end-t_cl)/del_t)+ceil(t_cl/del_t_fault);
x_sim=zeros(num_var,num_sim_idx); t_sim=zeros(1,num_sim_idx);
x_sim(:,1)=x_eq_pre; t_sim(1)=0; t_idx=1;

while t_sim(t_idx)<=t_end
    if t_sim(t_idx)<t_fault
        f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_pre_kron(x_sim(:,t_idx))+f_pre_kron(x));
        J=@(x) -eye(num_var)+del_t/2*J_pre_kron(x);
        t_sim(t_idx+1)=t_sim(t_idx)+del_t;
    elseif  t_sim(t_idx)<t_fault+t_cl
        f=@(x) -x+x_sim(:,t_idx)+del_t_fault/2*(f_fault_kron(x_sim(:,t_idx))+f_fault_kron(x));
        J=@(x) -eye(num_var)+del_t_fault/2*J_fault_kron(x);
        t_sim(t_idx+1)=t_sim(t_idx)+del_t_fault;
    else
        f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_post_kron(x_sim(:,t_idx))+f_post_kron(x));
        J=@(x) -eye(num_var)+del_t/2*J_post_kron(x);
        t_sim(t_idx+1)=t_sim(t_idx)+del_t;
    end
    x_sim(:,t_idx+1)=NR(f,J,x_sim(:,t_idx));
    t_idx=t_idx+1;
end
t_sim=t_sim(1:t_idx); x_sim=x_sim(:,1:t_idx);
tcomp_TDS=toc;
result_TDS=(max(x_sim(idx_delta,t_idx-1))-min(x_sim(idx_delta,t_idx-1)))<2*pi;