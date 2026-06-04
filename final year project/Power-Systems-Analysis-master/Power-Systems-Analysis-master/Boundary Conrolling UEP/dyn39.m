%% ================= Bus.con =================
% Column meanings:
% 1: Bus number
% 2: Voltage base (kV)
% 3: Bus type (1=SWING, 2=PV, 3=PQ)
% 4: Area number
% 5: Zone number
% 6: Status (1=in service, 0=out of service)
Bus.con = [ ... 
  1  16.5  1  0  4  1;
  2  18  1  0  5  1;
  3  13.8  1  0  3  1;
  4  230  1  0  2  1;
  5  230  1  0  2  1;
  6  230  1  0  2  1;
  7  230  1  0  2  1;
  8  230  1  0  2  1;
  9  230  1  0  2  1;
 ];

%% ================= Line.con =================
% Column meanings:
% 1: From bus
% 2: To bus
% 3: Base MVA
% 4: Voltage base (kV)
% 5: System frequency (Hz)
% 6: Status (1=in service)
% 7: Transformer tap ratio (0 or 1 if none)
% 8: Resistance R (p.u.)
% 9: Reactance X (p.u.)
% 10: Line charging B (p.u.)
% 11-15: Optional flags/padding
% 16: End-of-line / status
Line.con = [ ... 
  9  8  100  230  60  0  0  0.0119  0.1008  0.209  0  0  0  0  0  1;
  7  8  100  230  60  0  0  0.0085  0.072  0.149  0  0  0  0  0  1;
  9  6  100  230  60  0  0  0.039  0.17  0.358  0  0  0  0  0  1;
  7  5  100  230  60  0  0  0.032  0.161  0.306  0  0  0  0  0  1;
  5  4  100  230  60  0  0  0.01  0.085  0.176  0  0  0  0  0  1;
  6  4  100  230  60  0  0  0.017  0.092  0.158  0  0  0  0  0  1;
  2  7  100  18  60  0  0.07826087  0  0.0625  0  0  0  0  0  0  1;
  3  9  100  13.8  60  0  0.06  0  0.0586  0  0  0  0  0  0  1;
  1  4  100  16.5  60  0  0.07173913  0  0.0576  0  0  0  0  0  0  1;
 ];

%% ================= Breaker.con =================
% Column meanings:
% 1: Bus A
% 2: Bus B
% 3: Base MVA rating
% 4: Voltage rating (kV)
% 5: Frequency (Hz)
% 6: Status (1=closed)
% 7: Trip threshold multiplier
% 8: Thermal/current rating
% 9: Enable flag
% 10: Initial open/closed state
Breaker.con = [4  7  100  230  60  1  1.2  200  1  0;
                         4  5  100  230  60  1  1.2  200  1  0];

%% ================= Fault.con =================
% Column meanings:
% 1: Faulted bus
% 2: Base MVA
% 3: Voltage base (kV)
% 4: Frequency (Hz)
% 5: Fault type (1=3-phase)
% 6: Optional parameter
% 7: Fault resistance (p.u.)
% 8: Fault reactance (p.u.)
Fault.con = [
1   100   16.5   60  1  1.1  0  0.001 ;
2   100   18     60  1  1.1  0  0.001 ;
3   100   13.8   60  1  1.1  0  0.001 ;
4   100   230    60  1  1.1  0  0.001 ;
5   100   230    60  1  1.1  0  0.001 ;
6   100   230    60  1  1.1  0  0.001 ;
7   100  230  60  1  1.1  0  0.001];


%% ================= SW.con =================
% Slack/Swing Bus Data:
% 1: Bus number
% 2: Base MVA
% 3: Voltage base (kV)
% 4: Voltage magnitude setpoint
% 5: Voltage angle setpoint
% 6: Max active power
% 7: Min active power
% 8: Max voltage
% 9: Min voltage
% 10-13: Model/control flags
SW.con = [ ... 
  1  100  16.5  1.04  0  99  -99  1.1  0.9  0.8  1  1  1;
 ];

%% ================= PV.con =================
% Generator (PV) Data:
% 1: Bus number
% 2: Base MVA
% 3: Voltage base (kV)
% 4: Active power Pgen
% 5: Voltage setpoint
% 6: Max Q
% 7: Min Q
% 8: Max voltage
% 9: Min voltage
% 10-11: Control flags
PV.con = [ ... 
  2  100  18  1.63  1.025  99  -99  1.1  0.9  1  1;
  3  100  13.8  0.85  1.025  99  -99  1.1  0.9  1  1;
 ];

%% ================= PQ.con =================
% Load Data:
% 1: Bus number
% 2: Base MVA
% 3: Voltage base (kV)
% 4: Active power Pload
% 5: Reactive power Qload
% 6: Max voltage
% 7: Min voltage
% 8: Model flag
% 9: Status (1=in service)
PQ.con = [ ... 
  6  100  230  0.9  0.3  1.2  0.8  0  1;
  8  100  230  1  0.35  1.2  0.8  0  1;
  5  100  230  1.25  0.5  1.2  0.8  0  1;
 ];

%% ================= Syn.con =================
% Synchronous Machine Data:
% 1: Bus number
% 2: Base MVA
% 3: Voltage base (kV)
% 4: System frequency
% 5: Machine type / # poles
% 6-7: Optional parameters
% 8: d-axis synchronous reactance Xd
% 9: q-axis synchronous reactance Xq
% 10: d-axis transient reactance X'd
% 11: Subtransient reactance X''d
% 12-20: Time constants, damping, inertia H, other parameters
% 21+: Model flags for AVR/governor
Syn.con = [ ... 
  1  100  16.5  60  4  0  0  0.146  0.0608  0  8.96  0  0.0969  0.0969  0  0.31  0  47.28  0  0  0  1  1  0.002  0  0  1  1;
  2  100  18  60  4  0  0  0.8958  0.1198  0  6  0  0.8645  0.1969  0  0.535  0  12.8  0  0  0  1  1  0.002  0  0  1  1;
  3  100  13.8  60  4  0  0  1.3125  0.1813  0  5.89  0  1.2578  0.25  0  0.6  0  6.02  0  0  0  1  1  0.002  0  0  1  1;
 ];

%% ================= Exc.con =================
% Exciter Data:
% 1: Generator bus number
% 2-13: Exciter parameters (gain, time constants, limits, flags, saturation)
Exc.con = [ ... 
  1  2  5  -5  20  0.2  0.063  0.35  1  0.314  0.001  0.0039  1.555;
  2  2  5  -5  20  0.2  0.063  0.35  1  0.314  0.001  0.0039  1.555;
  3  2  5  -5  20  0.2  0.063  0.35  1  0.314  0.001  0.0039  1.555;
 ];

%% ================= Bus.names =================
% Bus names (strings for display / plots)
Bus.names = {... 
  'Bus 1'; 'Bus 2'; 'Bus 3'; 'Bus 4'; 'Bus 5'; 
  'Bus 6'; 'Bus 7'; 'Bus 8'; 'Bus 9'};
