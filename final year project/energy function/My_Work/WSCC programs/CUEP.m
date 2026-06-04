Y=[1.1411-2.298*1i 0.1323+0.7035*1i 0.1854+1.0611*1i
    0.1323+0.7035*1i 0.381-2.0202*1i 0.1965+1.2031*1i
    0.1854+1.0611*1i 0.1965+1.2031*1i 0.2723-2.3544*1i];
ths=[-0.1649 0.4987 0.2344];
M=[23.64 6.4 3.01];
for i=1:g
   for j=1:3
        C(i,j)=E(i)*E(j)*imag(Y(i,j));
        D(i,j)=E(i)*E(j)*real(Y(i,j)); 
   end
end
for k=1:g
   Pi(k)=Pm(k)-real(Y(k,k))*(abs(E(k)))^2; 
end
for it=1:1000
    Pcoi=0;
   for i=1:g
       Pcoi=Pcoi+Pi(i);
       for j=i+1:m
           Pcoi=Pcoi-2*D(i,j)*cos(th(i)-th(j));
       end
   end
   for i=1:g
       Pei(i)=0;
       for j=1:g
           if j~=i
               Pei(i)=Pei(i)+C(i,j)*sin(th(i)-th(j))+D(i,j)*cos(th(i)-th(j));
           end
       end
       f(i)=Pi(i)-Pei(i)-(H(i)/sum(H))*Pcoi;
   end
   fsum(it)=sum(f);
   th(it)=th;
end