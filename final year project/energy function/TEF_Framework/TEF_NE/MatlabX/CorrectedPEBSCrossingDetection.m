%PEBS crossing detection
%import Ybus
%Check the value of T, ths, sign of KE/PE
clc
tic
% ths=[-0.0199    0.2077    0.3073    0.2620    0.3868    0.2962    0.5595    0.2532 0.4642   -0.1625];
g=10;
% data=data1;

T=503; %simulation data length
del=data(1:T,1:g)*pi/180;
th=zeros(T,g);

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
Pi=Pm-(real(diag(Y1))').*((abs(E)).^2);

% fth=zeros(T,1);
% Pei=zeros(1,g);
% f=zeros(1,g);
% for t=1:T
%    Pcoi=0;
%    for i=1:g
%        Pcoi=Pcoi+Pi(i);
%        for j=i+1:g
%            Pcoi=Pcoi-2*D(i,j)*cos(th(t,i)-th(t,j));
%        end
%    end
%    for i=1:g
%        Pei(i)=0;
%        for j=1:g
%            if j~=i
%                Pei(i)=Pei(i)+C(i,j)*sin(th(t,i)-th(t,j))+D(i,j)*cos(th(t,i)-th(t,j));
%            end
%        end
%        f(i)=Pi(i)-Pei(i)-(H(i)/sum(H))*Pcoi;
%        fth(t)=fth(t)+((th(t,i)-ths(i))*f(i));
%    end
% end
% figure;plot(fth,'g');

% PEBS crossing point
% t=1;
% while t<500
%     if fth(t)>0
%         xx=t;
%         t=800;
%     end
%     t=t+1;
% end
% fth1=[fth(xx-1) fth(xx)];
% t1=[xx-1 xx];
% interp1(fth1,t1,0);
% dt=interp1(fth1,t1,0)-(xx-1);
% thpebs=th(xx-1,:)+(th(xx,:)-th(xx-1,:))*dt;
% %PEBS Crossing point end

Ws=314;
w=data(1:T,31:40)*Ws;
KE=zeros(T,1);
PE1=zeros(T,1);
PE2=zeros(T,1);
PE3=zeros(T,1);
PE4=zeros(T,1);
PE=zeros(T,1);
V=zeros(T,1);
% % w(100,:)=[0 0 0 0 0 0 0 0 0 0];
% % th(100,:)=[1.1132    1.7348    1.8764    1.7295    1.8472    1.7734    1.8268    1.5221    1.9174   -0.9541];
% % th(1,:)=CUEP;
% ths=CUEP;
wcoi(1:T,1)=(w(1:T,1)*H(1)+w(1:T,2)*H(2)+w(1:T,3)*H(3)+w(1:T,4)*H(4)+w(1:T,5)*H(5)+w(1:T,6)*H(6)+w(1:T,7)*H(7)+w(1:T,8)*H(8)+w(1:T,9)*H(9)+w(1:T,10)*H(10))/sum(H);
wc(1:T,1)=w(1:T,1)-wcoi;wc(1:T,2)=w(1:T,2)-wcoi;wc(1:T,3)=w(1:T,3)-wcoi;wc(1:T,4)=w(1:T,4)-wcoi;wc(1:T,5)=w(1:T,5)-wcoi;wc(1:T,6)=w(1:T,6)-wcoi;
wc(1:T,7)=w(1:T,7)-wcoi;wc(1:T,8)=w(1:T,8)-wcoi;wc(1:T,9)=w(1:T,9)-wcoi;wc(1:T,10)=w(1:T,10)-wcoi;
% tu=1402;% 
for t=1:T
    for i=1:g
        %KE(t)=KE(t)+0.5*(2*H(i)/Ws)*(wc(t,i))^2;
        PE1(t)=PE1(t)+Pi(i)*(th(t,i)-ths(i));
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
            PE2(t)=PE2(t)+C(i,j)*(cos(th(t,i)-th(t,j))-cos(ths(i)-ths(j)));
            PE3(t)=PE3(t)-D(i,j)*((th(t,i)-ths(i)+th(t,j)-ths(j))/(th(t,i)-th(t,j)-ths(i)+ths(j)))*(sin(th(t,i)-th(t,j))-sin(ths(i)-ths(j)));
            %PE3(t)=PE3(t)+D(i,j)*(ths(i)-th(t,i)+ths(j)-th(t,j))*cos((ths(i)-ths(j)+th(t,i)-th(t,j))/2)*(1-(((ths(i)-ths(j)-th(t,i)+th(t,j))/2)^2)/6+(((ths(i)-ths(j)-th(t,i)+th(t,j))/2)^4)/120);
            %PE4(t)=PE4(t)+sum(D(i,j)*(cos(th(t:tu,i)-th(t:tu,j))).*(th(t+1:tu+1,i)-th(t:tu,i)+th(t+1:tu+1,j)-th(t:tu,j)));
        end
    end
    PE(t)=(PE1(t)+PE2(t)+PE3(t));
    KE(t)=KE(t);
    V(t)=-KE(t)+PE(t); %for PEBS detection KE-PE, for calculating dV.. -KE+PE
end 
% hold on;plot(V,'r')
figure;plot(-KE)
hold on;plot(PE,'r')
%hold on;plot(PE4,'g')
% hold on;plot(PE,'g')
% V(100)+PE(120)
% KE(120)
% (V(1)+PE(2))/KE(2)

% %CCT estimation
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
% interp1(V1,t2,0)-2
dt=(interp1(V1,t2,0)-102)/100
toc