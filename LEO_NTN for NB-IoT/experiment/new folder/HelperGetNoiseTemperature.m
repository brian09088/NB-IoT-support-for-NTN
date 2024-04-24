function T = HelperGetNoiseTemperature(f,rx)
%HelperGetNoiseTemperature Obtain the receiver antenna noise temperature
%   T = HelperGetNoiseTemperature(FREQ,RX) returns the antenna noise
%   temperature of the receiver RX, corresponding to the frequency FREQ.
%
%   Example:
%   % Obtain the noise temperature of a receiver.
%
%   % Create a satelliteScenario object.
%   startTime = datetime(2024,4,18,22,0,0);
%   stopTime = startTime + days(1);
%   sampleTime = 60;
%   sc = satelliteScenario(startTime,stopTime,sampleTime);
%
%   % Add a satellite to the scenario.
%   semiMajorAxis = 6500000;             % m
%   eccentricity = 0;
%   inclination = 10;                     % degrees
%   rightAscensionOfAscendingNode = 0;    % degrees
%   argumentOfPeriapsis = 0;              % degrees
%   trueAnomaly = 0;                      % degrees
%   sat = satellite(sc,semiMajorAxis, ...
%       eccentricity, ...
%       inclination, ...
%       rightAscensionOfAscendingNode, ...
%       argumentOfPeriapsis, ...
%       trueAnomaly);
%
%   % Add a receiver to the scenario.
%   rx = receiver(sat);
%
%   % Obtain the noise temperature of the receiver antenna corresponding to
%   % a frequency of 30 GHz.
%   f = 30e9;
%   T = HelperGetNoiseTemperature(f,rx)

%   Copyright 2021 The MathWorks, Inc.

% Retrieve the receiver antenna object
rxAntenna = rx.Antenna;

% Calculate the maximum gain of the antenna
if isa(rxAntenna, 'satcom.satellitescenario.GaussianAntenna')
    g = pattern(rxAntenna, f, 0, 90); % dB
else
    gPattern = pattern(rxAntenna, f); % dB
    g = max(max(gPattern));
end


% Retrieve the receiver antenna gain to noise temperature ratio
% corresponding to the z-axis of the receiver
gByT = rx.GainToNoiseTemperatureRatio;

% Calculate the noise temperature in dBK
T_dB = g - gByT;  % dBK

% Calculate the noise temperature in K
T = 10^(T_dB/10); % K

end
