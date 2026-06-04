%Stable Equilibrium Point Kundur's Two Area System
%import Ybus,
tic
x0=data(100,1:10)*pi/180;
% x0=[-0.0140    0.2136    0.3139    0.2673    0.3936    0.3085    0.3547    0.2596    0.4706   -0.1553];
% E=[1.1113 1.1907 1.1777 1.1309 1.4476 1.2097 1.2809 1.1027 1.2569 1.0147];
 E=[1.08445585454352,1.11526151122246,1.11816916678907,1.06057082357884,1.24262796740079,1.16337714553885,1.16643634383566,1.05534775536221,1.09906529712836,1.03420238984330];
g=10;
C=zeros(g,g);
D=zeros(g,g);
for i=1:g
   for j=1:g
        C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j));
   end
end
% x0=th(172,:);
x0=MOD;
% x0=ths;
x0=x0-((x0(1)*H(1)+x0(2)*H(2)+x0(3)*H(3)+x0(4)*H(4)+x0(5)*H(5)+x0(6)*H(6)+x0(7)*H(7)+x0(8)*H(8)+x0(9)*H(9)+x0(10)*H(10))/sum(H))*ones(1,g);
A=[];
b=[];
Aeq=[];
beq=[];
lb=-pi*ones(1,g);
ub=pi*ones(1,g);
opts = optimset('Algorithm','interior-point');
% opts = optimset('Diplay','iter','Algorithm','sqp');
[x,fval,exitflag,~,lambda]=fmincon(@(x) 1,x0,A,b,Aeq,beq,lb,ub,@(x)SEPfunction(x,Pm,E,C,D,H,Y1),opts);
x=x-((x(1)*H(1)+x(2)*H(2)+x(3)*H(3)+x(4)*H(4)+x(5)*H(5)+x(6)*H(6)+x(7)*H(7)+x(8)*H(8)+x(9)*H(9)+x(10)*H(10))/sum(H))*ones(1,g);
ths=x
% ths=[-0.2167    0.4641    0.5568    0.6745    0.7452    0.5495    0.8665    0.6998    0.8306   -0.3024];
toc