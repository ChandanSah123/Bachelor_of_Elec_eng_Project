%Stable Equilibrium Point Kundur's Two Area System

%import Ybus,
clc;
H=[6.5 6.5 6.175 6.175];
E=[1.11337 1.111908 1.11164 1.10127];
Pm=[7.0008 7 7.19 7];
g=4;
C=zeros(g,g);
D=zeros(g,g);
for i=1:g
   for j=1:g
        C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
   end
end
x0=[1.3801    1.2140   -1.2768   -1.4539];
% x0=y(147,:);
A=[];
b=[];
Aeq=[];
beq=[];
lb=-pi*ones(1,g);
ub=pi*ones(1,g);
opts = optimset('Algorithm','interior-point');
% opts = optimset('Display','iter','Algorithm','sqp');
[x,fval,exitflag,~,lambda]=fmincon(@(x) 1,x0,A,b,Aeq,beq,lb,ub,@(x)SEPfunction(x,Pm,E,C,D,H,Y1),opts);
x=x-((x(1)*H(1)+x(2)*H(2)+x(3)*H(3)+x(4)*H(4))/sum(H))*ones(1,4)
x0