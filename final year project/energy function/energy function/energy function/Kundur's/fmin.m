%minimum ||f|| detection
%import Ybus
H=[6.5 6.5 6.175 6.175];%import inertia constant
T=800; %simulation data length
th=zeros(T,4);
Pm=[7.0008 7 7.19 7];
E=[1.11337 1.111908 1.11164 1.10127];
g=4;
C=zeros(g,g);
D=zeros(g,g);
Pi=[0 0 0 0];
for i=1:g
   for j=1:g
        C(i,j)=E(i)*E(j)*imag(Y(i,j));
        D(i,j)=E(i)*E(j)*real(Y(i,j)); 
   end
end
for k=1:g
   Pi(k)=Pm(k)-real(Y(k,k))*(abs(E(k)))^2; 
end
fth=zeros(T,1);
fsum=zeros(T,1);
f=[0 0 0 0];
Pei=[0 0 0 0];
for t=1:T
   th(t,:)=y(t,:);
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
       fsum(t)=fsum(t)+abs(f(i));
   end
end
figure;plot(fsum)
figure;plot(th)