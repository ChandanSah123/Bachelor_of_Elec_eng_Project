% Example 10.3   JHC June 5, 2017
% frequency deviation estimation 

t_more2 = [0:1/1440:(23+24*4)/1440];
X_mag = [];
X_phase = [];
for i = 1:24*4
  x = cos(2*pi*59*(t_more2(i:i+23)) + pi/4);
  X_phasor = fft(x);
  X_mag = [X_mag abs(X_phasor(2))];
  X_phase = [X_phase angle(X_phasor(2))*180/pi];
end
X_mag = X_mag/(24/2)/sqrt(2);

delta_f = (X_phase(3*24+1)-X_phase(1))*pi/180/(24*3/1440)/(2*pi)

% X_phase(3*24+1) = X_phase(73) = 24.3939, X_phase(1) = 42.5817

% delta_f = -1.0104  % not bad, exact answer would be -1 Hz

