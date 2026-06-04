%Stable Equilibrium Point Kundur's Two Area System
%import Ybus,
clc;
tic
x0=CUEP;
% x0=data(100,1:57)/57.32;
% x0=[-0.0185    0.2080    0.3090    0.2562    0.4005    0.3047    0.3518 0.2549    0.4654   -0.1529];
C=zeros(g,g);
D=zeros(g,g);
for i=1:g
   for j=1:g
        C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j));
   end
end
% x0=MOD;
x0=x0-((sum(x0.*H))/sum(H))*ones(1,g);
A=[];
b=[];
Aeq=[];
beq=[];
lb=-pi*ones(1,g);
ub=pi*ones(1,g);
opts = optimset('Algorithm','interior-point');
% opts = optimset('Diplay','iter','Algorithm','sqp');
[x,fval,exitflag,~,lambda]=fmincon(@(x) 1,x0,A,b,Aeq,beq,lb,ub,@(x)SEPfunction(x,Pm,E,C,D,H,Y1),opts);
% x=x-((x(1)*H(1)+x(2)*H(2)+x(3)*H(3)+x(4)*H(4)+x(5)*H(5)+x(6)*H(6)+x(7)*H(7)+x(8)*H(8)+x(9)*H(9)+x(10)*H(10))/sum(H))*ones(1,g);
CUEP=x
toc