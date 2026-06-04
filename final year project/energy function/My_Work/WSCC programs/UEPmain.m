clear all;
clc;
options = odeset('MaxStep',0.01,'InitialStep',0.01);
[T, Y]=ode45(@integrateth,[0 2], [-0.5916    1.1803    2.1368],options);
y=Y;