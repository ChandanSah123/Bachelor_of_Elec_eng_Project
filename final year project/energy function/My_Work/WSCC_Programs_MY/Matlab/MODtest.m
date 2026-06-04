% clear all;%MODTest.m
g=3;
ths=[-0.0651    0.2011    0.0837];
MOD=[2 3]';
[row col]=size(MOD);
H=[23.64 6.4 3.01];
H1=0;
H2=0;
th1=0;
th2=0;

for i=1:g
    c=0;
    for j=1:row
        if i==MOD(j)
            th1=th1+ths(i)*H(i);
            H1=H1+H(i);
            c=1;
        end
    end
    if c==0
        th2=th2+ths(i)*H(i);
        H2=H2+H(i);
    end
end

th1=th1/H1;
th2=th2/H2;
dth=th1-th2;
dth1=(pi-2*dth)*(H2/(H1+H2));
dth2=(pi-2*dth)*(H1/(H1+H2));

for i=1:g
    c=0;
    for j=1:row
        if i==MOD(j)
            thu(i)=ths(i)+dth1;
            c=1;
        end
    end
    if c==0
        thu(i)=ths(i)-dth2;
    end
end
MOD=thu