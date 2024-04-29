clc,clear,close all

startTime = datetime(2024,4,22,0,0,0);
stopTime = startTime + days(1);
sampleTime = 60;                                     % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);

%% tranfer by ST-2 TLE file
% ST-2
% 1 37606U 11022B   24100.53158358 -.00000219  00000+0  00000+0 0  9994
% 2 37606   0.0159 223.1612 0001592 155.1355  99.3470  1.00271446 47242

% Earth's gravitational parameter in m^3/s^2
mu = 3.986e14;

% Mean motion in radians per second
meanMotion = 1.00271446 * 2 * pi / 86400;

% Calculate the semi-major axis in meters
semiMajorAxis = (mu / (meanMotion^2))^(1/3);

% Orbital parameters from TLE
eccentricity = 0.0001592;
inclination = 0.0159;  % In degrees
rightAscensionOfAscendingNode = 223.1612;  % In degrees
argumentOfPeriapsis = 155.1355;  % In degrees
trueAnomaly = 99.3470;  % In degrees

% Create satellite in MATLAB
sat = satellite(sc, semiMajorAxis, eccentricity, inclination, ...
    rightAscensionOfAscendingNode, argumentOfPeriapsis, ...
    trueAnomaly, Name="ST-2");

gimbalrxSat = gimbal(sat);
gimbaltxSat = gimbal(sat);

gainToNoiseTemperatureRatio = 5;                                                        % dB/K
systemLoss = 3;                                                                         % dB
rxSat = receiver(gimbalrxSat,Name="Satellite Receiver",GainToNoiseTemperatureRatio= ...
    gainToNoiseTemperatureRatio,SystemLoss=systemLoss);

% s-band 2GHz
% txpower = 200mW = 23 dBm = -7 dBW
% bitRate ranging from 50 Mbps to 150 Mbps
frequency = 2e9;                                                                     % Hz
power = -7;                                                                            % dBW
bitRate = 100;                                                                         % Mbps
% systemLoss = PL;                                                                       % dB                                                                     % dB
systemLoss = 3;
txSat = transmitter(gimbaltxSat,Name="Satellite Transmitter",Frequency=frequency, ...
    power=power,BitRate=bitRate,SystemLoss=systemLoss);

% 594mm(sat Dish diameter)
dishDiameter = 5;                                                                    % meters
apertureEfficiency = 0.5;
gaussianAntenna(txSat,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);
gaussianAntenna(rxSat,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);

gs1 = groundStation(sc, 25.034, 121.564, Name="GS 1");
gs2 = groundStation(sc, 22.631, 120.302, Name="GS 2");

pointAt(gimbaltxSat,gs2);
pointAt(gimbalrxSat,gs1);

gimbalgs1 = gimbal(gs1);
gimbalgs2 = gimbal(gs2);

frequency = 2e9;                                                                          % Hz
power = -7;                                                                                % dBW
bitRate = 100;                                                                              % Mbps
txGs1 = transmitter(gimbalgs1,Name="GS 1 Transmitter",Frequency=frequency, ...
        Power=power,BitRate=bitRate);

requiredEbNo = 14;                                                                     % dB
rxGs2 = receiver(gimbalgs2,Name="Gs 2 Receiver",RequiredEbNo=requiredEbNo);

dishDiameter = 5;                                % meters
gaussianAntenna(txGs1,DishDiameter=dishDiameter);
gaussianAntenna(rxGs2,DishDiameter=dishDiameter);

pointAt(gimbalgs1,sat);
pointAt(gimbalgs2,sat);

lnk = link(txGs1,rxSat,txSat,rxGs2);
play(sc)
linkIntervals(lnk)