clc,clear,close all

startTime = datetime(2024,4,22,0,0,0);
stopTime = startTime + days(1);
sampleTime = 60;                                     % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);

%% tranfer by AsiaSAT-9 TLE file
% AsiaSAT-9
% 1 42942U 17057A   24102.05049493 -.00000367  00000-0  00000+0 0  9990
% 2 42942   0.0180 234.1068 0002305 140.3904 325.5366  1.00267846 23949

% Calculate semi-major axis from mean motion
mu = 3.986e14;  % Earth's gravitational parameter in m^3/s^2
n = 1.00267846 * 2 * pi / 86400;  % Mean motion in rad/s

semiMajorAxis = (mu / (n^2))^(1/3);  % Semi-major axis in meters
eccentricity = 0.0002305;
inclination = 0.0180;  % In degrees
rightAscensionOfAscendingNode = 234.1068;  % In degrees
argumentOfPeriapsis = 140.3904;  % In degrees
trueAnomaly = 325.5366;  % In degrees

% Create satellite in MATLAB
sat = satellite(sc, semiMajorAxis, eccentricity, inclination, ...
    rightAscensionOfAscendingNode, argumentOfPeriapsis, ...
    trueAnomaly, Name="AsiaSAT-9");

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