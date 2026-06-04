%% Settings
sys_case=9; % IEEE system
fault_line=4; % Line number of the faulted line
fault_frto_bus=1; % From or To bus where the fault is applied
t_cl=0.15; % clearing time

% Time domain simulation setting (TDS)
t_end=20; % end time for time domain simulation
del_t=0.05; % times tep for tds
del_t_fault=0.001; % times step for fault tds
del_t_BCU=0.01; % times step for BCU tds
t_fault=1; % time that the fault is applied

% Plot setting for phase portrait
plot_bus=2;