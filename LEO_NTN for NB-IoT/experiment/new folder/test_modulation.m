%% user defined data - Link Budget
Tant.Up = 290; % Antenna noise temperature in Kelvin for uplink , REMEMBER!! This is the satellite antenna
Tant.Down = 290; % Antenna noise temperature in Kelvin for downlink, REMEMBER!! This is the ground station antenna
Tant.Inter = 10; % Antenna noise temperature in Kelvin for interlink
NF = 2; % Noise figure in dB
T0 = 280; % Reference temperature in K
Pmin = -134; % Receiver sensitivity in dBW
freq.Inter = 2.2e9; % Frequency of interlink in Hz
freq.UpDown = 2.2e9; % Frequency of uplink/downlink in Hz
P_TR.Sat = 1; % Power in watts of the satellites
P_TR.GS = 5; % Power in watts of the ground stations
Bandwidth = 750e3; % Bandwidth in Hz
BitRate = 512e3; % Bit Rate in Hz (bit/s)

%% Uplink Modulation

% Modulation format : SC-FDMA
% Modulation order : QPSK/OPQSK

mod_order = 'OQPSK';

% SNR method
SNR_method = "berspec";

% BER specification
mod = mod_order;
% Bit Error Rate specification
ber_spec = 1e-6;

% Offset-Quadrature-Phase-Shift-Keying (OQPSK) is a variant of the QPSK modulation scheme 
% where the phase or timing of either the in-phase or 
% Quadrature component is shifted relative to each other by a one bit-period or 
% half a symbol-period Ts
Modulation(mod_order , ber_spec)

EbNo_spec = Modulation(mod, ber_spec) ;
SNR_spec = EbNo_To_SNR(EbNo_spec, BitRate, Bandwidth) ;


% This function computes the BER for different types of mod_orders. 
% For a given mod_order, mod_order order, and BER specification, it gives the required Eb/NO.

function [EbNo_spec] = Modulation(mod_order , ber_spec)

EbNo = 0:0.01:25;

switch mod_order
                    
    case {"OQPSK" ,"oqpsk"}
        dataenc = input('Differential encode (diff) or non differential (nondiff): ','s');
        BER = berawgn(EbNo, 'oqpsk' ,dataenc) ;

end

EbNo_spec = interp1(BER,EbNo, ber_spec) ;
semilogy (EbNo, BER)
hold on
semilogy (EbNo_spec,ber_spec, 'bp')
legend('BER', 'BER specification')
title ('BER vs Eb/No')
xlabel('Eb/No')
ylabel('BER')

end
