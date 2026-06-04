% clear all;
clc;
options = odeset('MaxStep',0.01,'InitialStep',0.01);
[T, Y]=ode45(@integratetheta,[0 2], [-1.0010   -1.0932    2.1478    0.0567],options);
y=Y;