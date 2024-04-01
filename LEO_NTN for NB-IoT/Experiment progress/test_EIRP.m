% EIRP calculation 公式驗算 & 實驗數據驗證

%% LEO Satellite parameters
% EIRP = 34 [dBW/MHz]

% Transmit Power: Given as 200 mW
% UE Antenna Gain: Given as 0 dBi
% Satellite Antenna Gain: Given as 30 dBi
% Bandwidth: Given as 2 GHz = 2000 Hz

%% Link budget simulation parameters
% Transmit power −6.99 dBW (200 mW -> 10log(0.2))
% Transmit antenna gain 0 dBi Antenna gain 
% EIRP (effective isotropic radiated power) −6.99 dBW (Transmit power + antenna_gain)

% antenna_gain = UE_antenna_gain + Sat_antenna_gain
EIRP = Transmit_Power + UE_antenna_gain + Sat_antenna_gain;
EIRP = 10 * log(200) + 0 + 30;

% EIRP per MHz = EIRP / Bandwidth(MHz) = 53 / 2000 = 0.0265 dBm/MHz

EIRP = 0.0265 - 30; 
EIRP = -29.9735;

% Finally, to express this in dBW/MHz, we subtract 30 dB
% dBm/MHz => dBW/MHz

% (~=34 dBW/MHz)