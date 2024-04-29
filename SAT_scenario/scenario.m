% add link analysis to transmitter

clc,clear,close all

startTime = datetime(2020,11,25,0,0,0);
stopTime = startTime + days(1);
sampleTime = 60;                                     % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);

semiMajorAxis = 10000000;                                                                  % meters
eccentricity = 0;
inclination = 60;                                                                          % degrees
rightAscensionOfAscendingNode = 0;                                                         % degrees
argumentOfPeriapsis = 0;                                                                   % degrees
trueAnomaly = 0;                                                                           % degrees
sat = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
        argumentOfPeriapsis,trueAnomaly,Name="Satellite");

gimbalrxSat = gimbal(sat);
gimbaltxSat = gimbal(sat);

gainToNoiseTemperatureRatio = 5;                                                        % dB/K
systemLoss = 3;                                                                         % dB
rxSat = receiver(gimbalrxSat,Name="Satellite Receiver",GainToNoiseTemperatureRatio= ...
    gainToNoiseTemperatureRatio,SystemLoss=systemLoss);

frequency = 27e9;                                                                     % Hz
power = 20;                                                                           % dBW
bitRate = 20;                                                                         % Mbps
systemLoss = 3;                                                                       % dB
txSat = transmitter(gimbaltxSat,Name="Satellite Transmitter",Frequency=frequency, ...
    power=power,BitRate=bitRate,SystemLoss=systemLoss);

dishDiameter = 0.5;                                                                    % meters
apertureEfficiency = 0.5;
gaussianAntenna(txSat,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);
gaussianAntenna(rxSat,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);

gs0 = groundStation(sc, 22.631, 120.302, Name="GS 0");
gs1 = groundStation(sc,Name="Ground Station 1");
latitude = 52.2294963;                                              % degrees
longitude = 0.1487094;                                              % degrees
gs2 = groundStation(sc,latitude,longitude,Name="Ground Station 2");

pointAt(gimbaltxSat,gs2);
pointAt(gimbalrxSat,gs1);

gimbalgs1 = gimbal(gs1);
gimbalgs2 = gimbal(gs2);

frequency = 30e9;                                                                          % Hz
power = 40;                                                                                % dBW
bitRate = 20;                                                                              % Mbps
txGs1 = transmitter(gimbalgs1,Name="Ground Station 1 Transmitter",Frequency=frequency, ...
        Power=power,BitRate=bitRate);

requiredEbNo = 14;                                                                     % dB
rxGs2 = receiver(gimbalgs2,Name="Ground Station 2 Receiver",RequiredEbNo=requiredEbNo);

dishDiameter = 5;                                % meters
gaussianAntenna(txGs1,DishDiameter=dishDiameter);
gaussianAntenna(rxGs2,DishDiameter=dishDiameter);

pointAt(gimbalgs1,sat);
pointAt(gimbalgs2,sat);

lnk = link(txGs1,rxSat,txSat,rxGs2);

linkIntervals(lnk)

play(sc)