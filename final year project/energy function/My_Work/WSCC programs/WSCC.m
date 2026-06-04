% data=data(:,2:2:24);
 %data=data(:,73:84);
Y=[1.1386-2.2966*1i 0.1290+0.7064*1i 0.1824+1.0637*1i
    0.1290+0.7064*1i 0.3744-2.0151*1i 0.1921+1.2067*1i
    0.1824+1.0637*1i 0.1921+1.2067*1i 0.2691-2.3516*1i];
% Y=[0.8450 - 2.9880i   0.2870 + 1.5130i   0.2100 + 1.2260i
 %   0.2870 + 1.5130i   0.4200 - 2.7240i   0.2130 + 1.0880i
  %  0.2100 + 1.2260i   0.2130 + 1.0880i   0.2770 - 2.3680i]; %no outage
ths=[-0.1782 0.5309 0.2711];
M=[23.64 6.4 3.01];
del=data(:,1:3)*pi/180;
T = 700;
for t=1:T
    dcoi=[del(t,1)*M(1)+del(t,2)*M(2)+del(t,3)*M(3)]/(sum(M));
    th(t,1)=del(t,1)-dcoi;
    th(t,2)=del(t,2)-dcoi;
    th(t,3)=del(t,3)-dcoi;
end
% th(1,:)=[-0.7968 2.1257 1.7381];
Pm = data(100, 4:6);
E=[1.054 1.050 1.017];

C12=E(1)*E(2)*imag(Y(1,2));
C13=E(1)*E(3)*imag(Y(1,3));
C23=E(2)*E(3)*imag(Y(2,3));
D12=E(1)*E(2)*real(Y(1,2));
D13=E(1)*E(3)*real(Y(1,3));
D23=E(2)*E(3)*real(Y(2,3));
for t=1:T
    P(1)=Pm(1)-real(Y(1,1))*E(1)^2;
    P(2)=Pm(2)-real(Y(2,2))*E(2)^2;
    P(3)=Pm(3)-real(Y(3,3))*E(3)^2;
    Pcoi=sum(P(:))-2*(D12*cos(th(t,1)-th(t,2))+D13*cos(th(t,1)-th(t,3))+D23*cos(th(t,2)-th(t,3)));
    f(1)=P(1)-(C12*sin(th(t,1)-th(t,2))+D12*cos(th(t,1)-th(t,2))+C13*sin(th(t,1)-th(t,3))+D13*cos(th(t,1)-th(t,3)))-(M(1)/(sum(M(:))))*Pcoi;
    f(2)=P(2)-(C12*sin(th(t,2)-th(t,1))+D12*cos(th(t,2)-th(t,1))+C23*sin(th(t,2)-th(t,3))+D23*cos(th(t,2)-th(t,3)))-(M(2)/(sum(M(:))))*Pcoi;
    f(3)=P(3)-(C13*sin(th(t,3)-th(t,1))+D13*cos(th(t,3)-th(t,1))+C23*sin(th(t,3)-th(t,2))+D23*cos(th(t,3)-th(t,2)))-(M(3)/(sum(M(:))))*Pcoi;
    fth(t)=f(1)*(th(t,1)-ths(1))+f(2)*(th(t,2)-ths(2))+f(3)*(th(t,3)-ths(3));
end
%Calculation of KE
w=data(:,10:12)*377;
%w(100,:)=[0 0 0];
%th(100,:)=[-0.7348 1.9721 1.5778];
for t=1:T
    wcoi=(w(t,1)*M(1)+w(t,2)*M(2)+w(t,3)*M(3))/(sum(M));
    KE(t)=0.5*(2*M(1)/377)*(w(t,1)-wcoi)^2+0.5*(2*M(2)/377)*(w(t,2)-wcoi)^2+0.5*(2*M(3)/377)*(w(t,3)-wcoi)^2;
    PE1(t)=(Pm(1)-real(Y(1,1))*E(1)^2)*(th(t,1)-ths(1))+(Pm(2)-real(Y(2,2))*E(2)^2)*(th(t,2)-ths(2))+(Pm(3)-real(Y(3,3))*E(3)^2)*(th(t,3)-ths(3));
    PE2(t)=C12*(cos(th(t,1)-th(t,2))-cos(ths(1)-ths(2)))+C13*(cos(th(t,1)-th(t,3))-cos(ths(1)-ths(3)))+C23*(cos(th(t,2)-th(t,3))-cos(ths(2)-ths(3)));
    PE3(t)=D12*((th(t,1)-ths(1)+th(t,2)-ths(2))/((th(t,1)-th(t,2))-(ths(1)-ths(2)))*((sin(th(t,1)-th(t,2))-sin(ths(1)-ths(2)))))+D13*((th(t,1)-ths(1)+th(t,3)-ths(3))/((th(t,1)-th(t,3))-(ths(1)-ths(3)))*((sin(th(t,1)-th(t,3))-sin(ths(1)-ths(3)))))+D23*((th(t,2)-ths(2)+th(t,3)-ths(3))/((th(t,2)-th(t,3))-(ths(2)-ths(3)))*((sin(th(t,2)-th(t,3))-sin(ths(2)-ths(3)))));
    PE(t)=PE1(t)+PE2(t)-PE3(t);
    V(t)=-PE(t)+KE(t);
end
%plot(fth)
hold on;
plot(V,'r')
hold on;
plot(-PE,'b')
hold on;
%plot(KE)


