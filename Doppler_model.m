%% functions & models for doppler shift effect setting (Ch2)

% doppler compensation is for problem statement (Ch3)

% 1.2.2.1
% Fading Due to Time Dispersion: Frequency-Selective Fading Channel
% 1.2.2.2 
% Fading Due to Frequency Dispersion: Time-Selective Fading Channel

function [fd] =  Doppler_model(fc, v, alpha)
    
    % fm = maximum doppler shift
    % Bd = doppler spectrum bandwidth
    % Bd = 2 * fm

    % Tc = coherence time, Tc ~= 1/fm
    % Bs : bandwidth period of the transmit signal
    % Ts : symbol period of the transmit signal
    % Bc : coherence bandwidth
    % sigma_tao denote : RMS delay spread
    c = 3 * (10^8);      % 光速 c = 3*(10^8) m/s

    fd = fc * ((v * cos(alpha))/c); % fd(kHz)

end
