%% link budget parameters

% This is the satellite antenna
Tant.Up = 290; % Antenna noise temperature in Kelvin for uplink
% This is the ground station antenna
Tant.Down = 290; % Antenna noise temperature in Kelvin for downlink 
Tant.Inter = 10; % Antenna noise temperature in Kelvin for interlink

NF = 2; % Noise figure in dB
T0 = 280; % Reference temperature in K
Pmin = -134; % Receiver sensitivity in dBW
freq.Inter = 2.2e9; % Frequency of interlink in Hz
freq.UpDown = 2.2e9; % Frequency of uplink/downlink in Hz
P_TR.Sat = 1; % Power in watts of the satellites
P_TR.GS = 5; % Power in watts of the ground stations
B = 750e3; % Bandwidth in Hz
BitRate = 512e3; % Bit Rate in Hz (bit/s)