% Example 11.1 Part 2   JHC June 6, 2017

% Two Area System from Sauer and Pai - Problem 8.1 and Example 11.1
% all em model

disp('Two-area power system power flow - base case')
% bus data format
% bus: 
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%       bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu

bus = [ 1  1.03       0   7.00   0.00  0.00   0.00  0.00  0.00 2  2.0  -2.0  22.0   1.5  .5;
	    2  1.01       0   7.00   0.00  0.00   0.00  0.00  0.00 2  2.0  -2.0  22.0   1.5  .5;
	    11 1.03       0   0.00   0.00  0.00   0.00  0.00  0.00 1  2.0  -2.0  22.0   1.5  .5;
	    12 1.01       0   7.00   0.00  0.00   0.00  0.00  0.00 2  2.0  -2.0  22.0   1.5  .5;
       101 1.00       0   0.00   0.00  0.00   0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
       102 1.00       0   0.00   0.00  0.00   0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
       111 1.00       0   0.00   0.00  0.00   0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
       112 1.00       0   0.00   0.00  0.00   0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
	    3  1.00       0   0.00   0.00 10.59  -0.735 0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
	    13 1.00       0   0.00   0.00 16.75  -0.899 0.00  0.00 3  0.0   0.0  230.0  1.5  .5]
       
% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio, tap phase, tapmax, tapmin, tapsize

line = [...
1   101 0.001   0.012    0.00    0.0  0. 0.  0.  0.;
2   102 0.001   0.012    0.00    0.0  0. 0.  0.  0.;
3   13  0.022   0.22     0.33    0.0  0. 0.  0.  0.;
3   13  0.022   0.22     0.33    0.0  0. 0.  0.  0.;
3   13  0.022   0.22     0.33    0.0  0. 0.  0.  0.;
3   102 0.002   0.02     0.03    0.0  0. 0.  0.  0.;
3   102 0.002   0.02     0.03    0.0  0. 0.  0.  0.;
11  111 0.001   0.012    0.00    0.0  0. 0.  0.  0.;
12  112 0.001   0.012    0.00    0.0  0. 0.  0.  0.;
13  112 0.002   0.02     0.03    0.0  0. 0.  0.  0.;
13  112 0.002   0.02     0.03    0.0  0. 0.  0.  0.;
101 102 0.005   0.05     0.075   0.0  0. 0.  0.  0.;
101 102 0.005   0.05     0.075   0.0  0. 0.  0.  0.;
111 112 0.005   0.05     0.075   0.0  0. 0.  0.  0.;
111 112 0.005   0.05     0.075   0.0  0. 0.  0.  0.]; 

% run power flow

[bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,30,1,'y',1);
