% function for calculate SNR

% Link Budget analysis

% input 
% param ----- unit
% signal bandwidth ---------------- kHz
% EIRP ---------------------------- dBW
% G/T  ---------------------------- dB/K
% signal attenuation : 
%       FSPL ---------------------- dB
%       scintillation loss (SL) --- dB
%       atmospheric loss (AL) ----- dB
%       shadow fading (SF) -------- dB
%       Boltzmann's constant k

% output 
% param ----- unit 
% SNR ----- dB

function [SNR] = calculate_SNR(signal_BW, EIRP, G_T, FSPL, SL, AL, SF)

    SNR = EIRP + G_T - 10*log10(signal_BW) - 10*log10(k) - FSPL;

end 