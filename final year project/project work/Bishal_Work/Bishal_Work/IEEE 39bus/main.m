%% ------------------------------------------------------Main Program-------------------------------------------------%%
% open IEEE39 bus system
clc
clear all
open_system('D:\Final year project\Final_Year_Project_BEL[1]\Bishal_Work\IEEE 39bus\Simulink_model.slx');    % open 39bus system
[line_data, trans_data, load_data, generator_data1, generator_data2, mechanical_data] = extract_data();       % call extract_data function 
Ybus=Ybus_formation(line_data,trans_data);                                                                    % call Ybus formation function

gen_bus = [];
xd1=[];        % xd' = transient reactance for transient analysis
xd2=[];           % xd" = sub-transient reactant for short circuit analysis

for i = 1:length(generator_data2)
    gen_bus(i)=i;
    xd1(i)= generator_data2(i).Xdd;
    xd2(i)= generator_data2(i).Xddd;

end

Yextented_transient=build_extended_Ybus(Ybus,gen_bus,xd1);
Yextented_shortCircuit=build_extended_Ybus(Ybus,gen_bus,xd2);

%disp(Ybus);
