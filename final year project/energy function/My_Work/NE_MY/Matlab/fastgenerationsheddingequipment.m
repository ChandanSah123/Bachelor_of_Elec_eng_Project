g=zeros(10,2);
angle=zeros(10,2);
for g=1:10
    tt=[135 137 140];
    t2=0;t1=0.02;t0=0.05;t=0.2;
    w2=data(tt(1),30+g); w1=data(tt(2),30+g); w0=data(tt(3),30+g);
    del2=data(tt(1),g); del1=data(tt(2),g); del0=data(tt(3),g);
    alpha1=(w2-w1)/(t2-t1);
    alpha0=(w1-w0)/(t1-t0);
    alpha2=(alpha1-alpha0)/(t2-t0);
    w3=w2+alpha1*(t-t2)+alpha2*(t-t1)*(t-t2);
    del=del0+w2*(t-t0)+alpha1*((t^2-t0^2)/2-t2*(t-t0))+alpha2*((t^3-t0^3)/3-(t1+t2)*(t^2-t0^2)/2+t1*t2*(t-t0));
    w(g,1:2)=[data(tt(1)+20,30+g) w3];
    angle(g,1:2)=[data(tt(1)+20,g) del];
end
w
angle