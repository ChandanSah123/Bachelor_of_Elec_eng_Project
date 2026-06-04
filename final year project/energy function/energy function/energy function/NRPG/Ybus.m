%Formation of Admittance matrix, Kundur's Two Area System
%Update E, Pm, Xd', Load, line outage, shed, H as applicable
tic
clc;
global Y1
global shed
global remain
global Pm
global H
global E
global g
line=line1;
E=[1.0100    1.0406    1.0391    1.0628    1.0842    1.0476    1.0123    1.0806    1.0709    1.0137    1.0275    1.0002    1.0072    1.0036    0.9678    0.9879    0.9876    0.9784    0.9882    0.9838    0.9866    0.9882    0.9882    0.9887    0.9887    0.9862    1.0652    1.0532    1.0276    1.0273    1.0469    1.0328    1.0142    0.9924   0.9935    1.0537    1.0533    1.0597    0.9720    0.9896    0.9948    0.9216    1.0046    0.9876    1.0036    1.0032    1.0041    1.0025    1.0803    1.0494    1.0284    1.0305    1.0019    1.0002    1.0014    1.0156    0.9643];
for k=1:g
    H(k)=h(k)*MVA_base(k)/100;
end
Pm=data1(50,g+1:g*2);

shg=[31]; %shedded generator
shed=0.9;
remain=1-shed;
Pm(shg)=Pm(shg)*remain;
H(shg)=H(shg)*remain;
line(372+shg,4)=line(372+shg,4)*(1+shed/remain);

shg=[31]; %shedded generator
shed=0.0;
remain=1-shed;
Pm(shg)=Pm(shg)*remain;
H(shg)=H(shg)*remain;
line(372+shg,4)=line(372+shg,4)*(1+shed/remain);

bus=303;
Y1=zeros(bus,bus);
nload=246;
for k=1:429
   z=line(k,3)+1i*line(k,4);
   Y1(line(k,1),line(k,2))=-(1/z);
   Y1(line(k,2),line(k,1))=-(1/z);
   Y1(line(k,1),line(k,1))=Y1(line(k,1),line(k,1))+(1/z)+1i*line(k,5)*0.5;
   Y1(line(k,2),line(k,2))=Y1(line(k,2),line(k,2))+(1/z)+1i*line(k,5)*0.5;
end
%conversion of loads as admittance load
for k=1:nload
   Y1(load(k,1),load(k,1))=Y1(load(k,1),load(k,1))+load(k,3)/((load(k,2))^2)-1i*load(k,4)/((load(k,2))^2);
end
Ynn=Y1(247:303,247:303);
Ynr=Y1(247:303,1:246);
Yrn=Y1(1:246,247:303);
Yrr=Y1(1:246,1:246);
Yred=Ynn-((Ynr/Yrr)*Yrn);
Y1=Yred;

th=ang;
vn=zeros(1,g);
for k=1:g
   vn(k)=E(k)*(cosd(th(k))+1i*sind(th(k)));
end
Pe=zeros(1,g);
for k=1:g
    Pe(1,k)=(real(Yred(k,k))*E(k)^2);
    for n=1:g
        if k~=n
            Pe(1,k)=Pe(1,k)+abs(vn(k))*abs(vn(n))*abs(Yred(k,n))*cos(angle(vn(k))-angle(vn(n))-angle(Yred(k,n)));
        end
    end
end
Pm-Pe
toc