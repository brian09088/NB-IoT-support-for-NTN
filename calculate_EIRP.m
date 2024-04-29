%% functions for calculate EIRP

% Satellite effective isotropic radiated power (dB/MHz)
% 衛星有效各向同性輻射功率

%% LEO Satellite parameters
% EIRP = 34 [dBW/MHz]

% Transmit Power: Given as 200 mW
% UE Antenna Gain: Given as 0 dBi
% Satellite Antenna Gain: Given as 30 dBi
% Bandwidth: Given as 2 GHz = 2000 Hz

%% Link budget simulation parameters
% Transmit power = −6.99 dBW (200 mW -> 10log(0.2))
% Transmit antenna gain = 0 dBi (Antenna gain) 
% EIRP (effective isotropic radiated power) −6.99 dBW (Transmit power + antenna_gain)

function [EIRP] = calculate_EIRP(UE_antenna_power, UE_antenna_gain, Sat_antenna_gain)
    
    UE_antenna_gain = 0;
    Sat_antenna_gain = 0;
    UE_antenna_W = UE_antenna_power / 1000;
    UE_antenna_dBW = 10 * log10(UE_antenna_W);
    % Transmit power = 200mW ~= 23 dBm ~= −6.99 dBW
    % dBW = dBm -30
    
    EIRP = UE_antenna_dBW + UE_antenna_gain + Sat_antenna_gain; 

end
