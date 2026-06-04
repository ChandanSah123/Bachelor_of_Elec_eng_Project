% clear all;%MOD Test
tic
g=10;
% ths=[-0.2167    0.4641    0.5568    0.6745    0.7452    0.5495    0.8665    0.6998    0.8306   -0.3024];%2 axis model
ths=[-0.0199    0.2077    0.3073    0.2620    0.3868    0.2962    0.5595    0.2532 0.4642   -0.1625];
MOD1=[1:9]';
MOD=MOD1;
[row col]=size(MOD);
 H=[42 30.3 35.8 28.6 26 34.8 26.4 24.3 34.5 500];
% H(shg)=H(shg)*remain;
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
toc