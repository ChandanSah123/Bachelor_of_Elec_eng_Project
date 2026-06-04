%PEBS crossing detection
%import Ybus
clc
g=10;
T=5000; %simulation data length
del=data(1:T,1:g)*pi/180;
th=zeros(T,g);
Ws=314;
w=data(1:T,31:40)*Ws;
dcoi=(del(:,1)*H(1)+del(:,2)*H(2)+del(:,3)*H(3)+del(:,4)*H(4)+del(:,5)*H(5)+del(:,6)*H(6)+del(:,7)*H(7)+del(:,8)*H(8)+del(:,9)*H(9)+del(:,10)*H(10))/(sum(H));
    th(:,1)=del(:,1)-dcoi;
    th(:,2)=del(:,2)-dcoi;
    th(:,3)=del(:,3)-dcoi;
    th(:,4)=del(:,4)-dcoi;
    th(:,5)=del(:,5)-dcoi;
    th(:,6)=del(:,6)-dcoi;
    th(:,7)=del(:,7)-dcoi;
    th(:,8)=del(:,8)-dcoi;
    th(:,9)=del(:,9)-dcoi;
    th(:,10)=del(:,10)-dcoi;

Pi=zeros(1,g);
% Edd=[1.0934    1.1826    1.1496    1.0829    1.3974    1.1989    1.1825    1.0713    1.1373     1.0365]; %Absolute value, magnitude only
% Eqd=[1.0938    1.1828    1.1498    1.0830    1.3976    1.1991    1.1825    1.0713    1.1374     1.0368];
% alpha=zeros(g,g);
% beta=zeros(g,g);
% for i=1:g
%     for j=1:i
%         alpha(i,j)=Edd(i)*Edd(j)+Eqd(i)*Eqd(j);
%         beta(i,j)=Edd(i)*Eqd(j)-Eqd(i)*Edd(j);
%         alpha(j,i)=alpha(i,j);
%         beta(j,i)=-beta(i,j);
%     end
% end
% Pi=Pm'-diag(alpha).*diag(G);
KE=zeros(T,1);
PE1=zeros(T,1);
PE2=zeros(T,1);
PE3=zeros(T,1);
PE=zeros(T,1);
V=zeros(T,1);
wcoi(1:T,1)=(w(1:T,1)*H(1)+w(1:T,2)*H(2)+w(1:T,3)*H(3)+w(1:T,4)*H(4)+w(1:T,5)*H(5)+w(1:T,6)*H(6)+w(1:T,7)*H(7)+w(1:T,8)*H(8)+w(1:T,9)*H(9)+w(1:T,10)*H(10))/sum(H);
wc(1:T,1)=w(1:T,1)-wcoi;wc(1:T,2)=w(1:T,2)-wcoi;wc(1:T,3)=w(1:T,3)-wcoi;wc(1:T,4)=w(1:T,4)-wcoi;wc(1:T,5)=w(1:T,5)-wcoi;wc(1:T,6)=w(1:T,6)-wcoi;
wc(1:T,7)=w(1:T,7)-wcoi;wc(1:T,8)=w(1:T,8)-wcoi;wc(1:T,9)=w(1:T,9)-wcoi;wc(1:T,10)=w(1:T,10)-wcoi;
% ths=CUEP;
Eqdp=[1.1116    1.1909    1.1778    1.1310    1.4478    1.2099    1.2811    1.1028    1.2571    1.0151]; %emf at the peak point
Eddp=[1.1113    1.1907    1.1777    1.1309    1.4476    1.2097    1.2809    1.1027    1.2569    1.0147]; %emf at the peak point
Edd1(1000,:)=Eddp;
Eqd1(1000,:)=Eqdp;
for t=1:T
    Edd=(abs(Edd1(t,:))'+Eddp'); %Absolute value, magnitude only
    Eqd=(abs(Eqd1(t,:))'+Eqdp');
    alpha=zeros(g,g);
    beta=zeros(g,g);
    for i=1:g
        for j=1:i
        alpha(i,j)=Edd(i)*Edd(j)+Eqd(i)*Eqd(j);
        beta(i,j)=Edd(i)*Eqd(j)-Eqd(i)*Edd(j);
        alpha(j,i)=alpha(i,j);
        beta(j,i)=-beta(i,j);
        end
    end
    Pi=Pm'-diag(alpha).*diag(G);
    
    for i=1:g
        KE1(t,i)=0.5*(2*H(i)/Ws)*(wc(t,i))^2;
        PE1(t)=PE1(t)-Pi(i)*(th(t,i)-ths(i));
    end
    wcr=0;
    wsys=0;
    Hcr=0;
    Hsys=0;
    MOD=MOD1;
    for i=1:g
    c=0;
    for j=1:numel(MOD)
        if i==MOD(j)
            wcr=wcr+H(i)*wc(t,i);
            Hcr=Hcr+H(i);
            c=1;
        end
    end
    if c==0
        wsys=wsys+H(i)*wc(t,i);
        Hsys=Hsys+H(i);
    end
    end
    wcr=wcr/Hcr;
    wsys=wsys/Hsys;
    weq(t)=wcr-wsys;
    Heq=Hcr*Hsys/(Hcr+Hsys);
    KE(t)=0.5*(2*Heq/Ws)*(weq(t))^2;
    for i=1:g-1
        for j=i+1:g
            x1=-cos(th(t,i)-th(t,j))+cos(ths(i)-ths(j));
            x2=sin(th(t,i)-th(t,j))-sin(ths(i)-ths(j));
            x3=(th(t,i)-ths(i)+th(t,j)-ths(j))/(th(t,i)-th(t,j)-ths(i)+ths(j));
            x4=(th(t,i)-ths(i)+th(t,j)-ths(j))/(th(t,i)-th(t,j)-ths(i)+ths(j));
            
            PE2(t)=PE2(t)+beta(i,j)*alpha(i,j)*(x1)-beta(i,j)*B(i,j)*(x2);
            PE3(t)=PE3(t)+G(i,j)*alpha(i,j)*(x3)*(x2)+G(i,j)*beta(i,j)*x4*(-x1);
        end
    end
    PE(t)=(PE1(t)+PE2(t)+PE3(t));
    KE(t)=KE(t);
    V(t)=-KE(t)+PE(t); %for PEBS detection KE-PE, for calculating dV.. -KE+PE
end
hold on;plot(PE,'r')
t=1;
while t<T
    if V(t)<0
        yy=t;
        t=T+10;
    end
    t=t+1;
end
V1=[V(yy-1) V(yy)];
t2=[yy-1 yy];
dt=(interp1(V1,t2,0)-1003)/1000