%PEBSCrossingDetection.m
%import Ybus
g=3;
%ths=[-0.0651    0.2011    0.0837];
 %ths =[0.0077    0.0100   -0.0817];
 ths =[-0.1782    0.5309    0.2711];
H=[23.64 6.4 3.01];
%T=length(data(:,1));%simulation data length
T=300;
th=zeros(T,g);
%data=data1;
Y1=Yint_post;
Pm=data(100,4:6);
E=[1.0728 1.0775 1.0609];
del=data(:,1:3);
for t=1:T
    dcoi=(del(t,1)*H(1)+del(t,2)*H(2)+del(t,3)*H(3))/(sum(H));
    th(t,1)=del(t,1)-dcoi;
    th(t,2)=del(t,2)-dcoi;
    th(t,3)=del(t,3)-dcoi;
end
C=zeros(g,g);
D=zeros(g,g);
for i=1:g
   for j=1:g
        C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
   end
end
Pi=zeros(1,g);
for k=1:g
   Pi(k)=Pm(k)-real(Y1(k,k))*(abs(E(k)))^2; 
end
fth=zeros(T,1);
Pei=zeros(1,g);
f=zeros(1,g);
for t=1:T
   Pcoi=0;
   for i=1:g
       Pcoi=Pcoi+Pi(i);
       for j=i+1:g
           Pcoi=Pcoi-2*D(i,j)*cos(th(t,i)-th(t,j));
       end
   end
   for i=1:g
       Pei(i)=0;
       for j=1:g
           if j~=i
               Pei(i)=Pei(i)+C(i,j)*sin(th(t,i)-th(t,j))+D(i,j)*cos(th(t,i)-th(t,j));
           end
       end
       f(i)=Pi(i)-Pei(i)-(H(i)/sum(H))*Pcoi;
       fth(t)=fth(t)+((th(t,i)-ths(i))*f(i));
   end
end
figure;plot(fth,'g');
w=data(:,10:12)*377;
KE=zeros(T,1);
PE1=zeros(T,1);
PE2=zeros(T,1);
PE3=zeros(T,1);
PE=zeros(T,1);
V=zeros(T,1);
w(100,:)=[0 0 0];
%th(100,:)=[-0.7964    2.1263    1.7339];

for t=1:T
    wcoi=(w(t,1)*H(1)+w(t,2)*H(2)+w(t,3)*H(3))/(sum(H));
    for i=1:g
%         KE(t)=KE(t)+0.5*(2*H(i)/377)*(w(t,i)-wcoi)^2;
        wc(i)=(w(t,i)-wcoi);
        PE1(t)=PE1(t)+Pi(i)*(th(t,i)-ths(i));
    end
    Hcr=H(3)+H(2);
    Hsys=H(1);
    wcr=(H(2)*wc(2)+H(3)*wc(2))/Hcr;
    wsys=(H(1)*wc(1))/Hsys;
    weq=-wcr+wsys;
    Heq=Hcr*Hsys/(Hcr+Hsys);
    KE(t)=0.5*(2*Heq/377)*(weq)^2;
    for i=1:g-1
        for j=i+1:g
            PE2(t)=PE2(t)+C(i,j)*(cos(th(t,i)-th(t,j))-cos(ths(i)-ths(j)));
            PE3(t)=PE3(t)-D(i,j)*((th(t,i)-ths(i)+th(t,j)-ths(j))/(th(t,i)-th(t,j)-ths(i)+ths(j)))*(sin(th(t,i)-th(t,j))-sin(ths(i)-ths(j)));
        end
    end
    PE(t)=PE1(t)+PE2(t)+PE3(t);
    V(t)=KE(t)-PE(t); %for PEBS detection KE-PE, for calculating dV.. -KE+PE
end
figure;plot(V,'r')
hold on;plot(-PE,'b')