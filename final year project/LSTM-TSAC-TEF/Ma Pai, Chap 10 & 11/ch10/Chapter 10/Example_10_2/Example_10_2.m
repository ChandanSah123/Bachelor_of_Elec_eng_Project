% Example 10.2: 3-pt filter 
jay = sqrt(-1);

Xmag = [0.7073; 0.7075; 0.7065];
Xph  = [59.7593; 119.5788; 179.4992]*pi/180;
Xph_new = [59.7593; 119.5788-60; 179.4992-120]*pi/180;

X0p = sum(Xmag.*exp(jay*Xph))/2
abs(X0p), angle(X0p)*180/pi

X0pmag = sum(Xmag)/3
X0pph  = sum(Xph)*180/pi/3 - 60

X0p_new = sum(Xmag.*exp(jay*Xph_new))/3
X0p_newmag = abs(X0p_new)
X0p_newph = angle(X0p_new)*180/pi

Pmag = sin(24/2*(-0.1)/1440*2*pi)/(24*sin(1/2*(-0.1)/1440*2*pi))

Pph = (24-1)*(-0.1)/2/1440*2*pi

Pph = (24-1)*(-0.1)/2/1440*360

Xph_est = 59.7593 - Pph

% Pph = (24-1)*(-0.1)/2/1440  % Dotta's value

% Pph = (24-1)*(-0.1)/2/1440*2*pi/180 

% Results
% X0p =  -0.3497 + 0.6163i
% ans =    0.7086
% ans =  119.5761
% X0pmag =    0.7071
% X0pph =   59.6124
% X0p_new =   0.3577 + 0.6100i
% X0p_newmag =    0.7071
% X0p_newph =   59.6125
% Pmag =    1.0000
% Pph =   -0.0050
% Pph =   -0.2875
% Xph_est =   60.0468
