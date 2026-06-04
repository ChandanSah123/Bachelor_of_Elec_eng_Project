% Plot file for Example 10.4
% The radial chain of multiple-machine system to illustrate 
% electromechanical wave propagation 
% base MVA 100
% data_emwave_1.m

% use PST function s_simu to simulate the 12-machine system in the data file 
% data_emwave_1.m


% commands for plotting
set_font_size

% one-second plot
figure, plot(t(1:141),mac_spd(1:11,1:141),[0 1],[0.9999 0.9999])
xlabel('Time (sec)'), ylabel('Frequency(pu)')

% five-second plot
figure, plot(t,mac_spd(1:11,:))
xlabel('Time (sec)'), ylabel('Frequency(pu)')

reset_font_size

figure, mesh(t,[1:1:11],mac_spd(1:11,:)) % use rotate button for a better display 
                                  % of the time response 
xlabel('Time (sec)'), ylabel('machine number'), zlabel('Machine Speed (pu)')




