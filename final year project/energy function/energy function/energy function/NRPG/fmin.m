tic
E=[1.0840    1.1150    1.1174    1.0593    1.2412    1.1571    1.1294    1.0549    1.0985    1.0342];%original
%H=[42 30.3 35.8 28.6 26 34.8 26.4 24.3 34.5 500];
g=10; %no of generators
%Pm=[250.0 520.8 650.0 632.0 508.0 650.0 560.0 540.0 830.0 1000.0]/100;
T1=200;
Pi=zeros(g,1);
C=zeros(g,g);
D=zeros(g,g);
for i=1:g
    for j=1:g
        C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
    end
end
for k=1:g
    Pi(k)=Pm(k)-real(Y1(k,k))*(abs(E(k)))^2;
end
fsum=zeros(T1,1);
f=zeros(g,1);
Pei=zeros(g,1);
for it=1:length(y)
Pcoi=0;
fsum(it)=0;
th=y(it,:)';
th=th-ones(g,1)*((th(1)*H(1)+th(2)*H(2)+th(3)*H(3)+th(4)*H(4)+th(5)*H(5)+th(6)*H(6)+th(7)*H(7)+th(8)*H(8)+th(9)*H(9)+th(10)*H(10))/sum(H));
th=th';
    for i=1:g
       Pcoi=Pcoi+Pi(i);
       for j=i+1:g
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
        fsum(it)=fsum(it)+abs(f(i));
    end
end
figure;plot(fsum)
toc