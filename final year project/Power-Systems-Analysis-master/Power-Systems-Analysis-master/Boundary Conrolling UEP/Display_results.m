%% Display results
figure;
subplot(2,1,1); hold all; grid on; box on;
plot(t_sim,x_sim(idx_delta,:))
xlim([0 t_end]); set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('\delta');
subplot(2,1,2); hold all; grid on; box on;
plot(t_sim,x_sim(idx_omega,:))
xlim([0 t_end]); set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('\omega');

figure;
subplot(3,1,1); hold all; grid on; box on;
plot(t_fault_sim,V_ke,'r--');
plot(t_fault_sim,V_p+V_d,'b--');
plot(t_fault_sim,V_total);
% xlim([t_fault t_fault+t_cl])
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('Energy function');
legend('V_{ke}','V_{pe}','V_{ke}+V_{pe}')

subplot(3,1,2); hold all; grid on; box on;
plot(t_fault_sim,PEBS);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)');  ylabel('PEBS function');

subplot(3,1,3); hold all; grid on; box on;
plot(t_bcu_sim,theta_abs);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('iteration');  ylabel('|f|');

% Phase Portrait
figure;
plot_axis=[plot_bus plot_bus+num_gen];
plot_rng=[-pi pi -10 10];
x_plot_pre=x_eq_pre;
x_plot_post=x_eq_post;
[x1_plot_ij,x2_plot_ij]=meshgrid(linspace(plot_rng(1),plot_rng(2),30),linspace(plot_rng(3),plot_rng(4),30));
for i=1:size(x1_plot_ij,1)
    for j=1:size(x2_plot_ij,2)
        x_plot_pre(plot_axis)=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        x_plot_post(plot_axis)=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        dx_plot_pre=f_fault_kron(x_plot_pre);
        dx_plot_post=f_post_kron(x_plot_post);
        dx1_plot_pre(i,j)=dx_plot_pre(plot_axis(1));
        dx2_plot_pre(i,j)=dx_plot_pre(plot_axis(2));
        dx1_plot_post(i,j)=dx_plot_post(plot_axis(1));
        dx2_plot_post(i,j)=dx_plot_post(plot_axis(2));
        num_prevent_zero=E_kron*(x_plot_post(idx_delta)-x_eq_post(idx_delta));
        num_prevent_zero(num_prevent_zero==0)=inf;
        V_plot(i,j)=0.5*sum(M_gen.*x_plot_post(idx_omega).^2)-P_post_kron'*(x_plot_post(idx_delta)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_plot_post(idx_delta))-cos(E_kron*x_eq_post(idx_delta)))+((abs(E_kron)*(x_plot_post(idx_delta)-x_eq_post(idx_delta)))./(num_prevent_zero))'*diag(gij_post_kron)*(sin(E_kron*x_plot_post(idx_delta))-sin(E_kron*x_eq_post(idx_delta)));
    end
end
subplot(1,2,1); hold all; grid on; box on;
plot(x_sim(plot_axis(1),idx_tf:idx_tcl),x_sim(plot_axis(2),idx_tf:idx_tcl),'LineWidth',3);
streamslice(x1_plot_ij,x2_plot_ij,dx1_plot_pre,dx2_plot_pre);
axis(plot_rng); 
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); title('Fault-on Dynamics')
subplot(1,2,2); hold all; grid on; box on;
scatter(theta_u(plot_axis(1)),0,100,'r','filled')
scatter(x_fault(plot_axis(1),PEBS_idx),0,100,'b','filled')
plot(x_sim(plot_axis(1),:),x_sim(plot_axis(2),:),'b','LineWidth',2);
streamslice(x1_plot_ij,x2_plot_ij,dx1_plot_post,dx2_plot_post);
axis(plot_rng); 
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); title('Post-contingency Dynamics'); legend('u.e.p. (BCU)','u.e.p. (PEBS)')

% Plot Energy function
figure; hold all; grid on; box on;
scatter(theta_u(plot_axis(1)),0,100,'r','filled')
scatter(x_fault(plot_axis(1),PEBS_idx),0,100,'b','filled')
streamslice(x1_plot_ij,x2_plot_ij,dx1_plot_post,dx2_plot_post);
contour(x1_plot_ij,x2_plot_ij,V_plot,'LineWidth',2);
mesh(x1_plot_ij,x2_plot_ij,V_plot);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); zlabel('V'); legend('u.e.p. (BCU)','u.e.p. (PEBS)')

% Display results
disp(['System: ' num2str(sys_case) ' bus system with fault on line ' num2str(fault_line) ' and bus ' num2str(line_frto(fault_line,fault_frto_bus))])
disp(['Result:                        (TDS) ' 'un'*~result_TDS 'stable  (BCU) ' 'un'*~result_BCU  'stable  (PEBS) ' 'un'*~result_PEBS 'stable'])
disp(['V critical:                    (@t_cl) ' num2str(V_total(idx_tcl),'%.3f') '     (BCU) ' num2str(Vcr_BCU,'%.3f')  '     (PEBS) ' num2str(Vcr_PEBS,'%.3f')])
disp(['Critical Clearing Time: (t_cl) ' num2str(t_cl,'%.3f') '       (BCU) ' num2str(tcr_BCU,'%.3f') '     (PEBS) ' num2str(tcr_PEBS,'%.3f')])
disp(['Computation Time:      (TDS) ' num2str(tcomp_TDS,'%.2f') '       (BCU) ' num2str(tcomp_BCU,'%.2f')  '       (PEBS) ' num2str(tcomp_PEBS,'%.2f')])

