% clear all;%MOD Test
g=4;
ths=[0.3165    0.1504   -0.1572   -0.3343]; %post fault SEP
MOD=[1 2]'; %give MOD generator numbers
[row col]=size(MOD);
H=[6.5 6.5 6.175 6.175]; %give H values
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