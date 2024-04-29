clc,clear,close all
%% Interference from Satellite Constellation on Communications Link

% Create Satellite Scenario
startTime = datetime(2024,4,18,22,0,0);              
stopTime = startTime + minutes(10);                    
sampleTime = 60;                             % In s
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Add Medium-Earth Orbit Satellite
semiMajorAxis = 12000000;                    % In m
eccentricity = 0;
inclination = 8;                             % In degrees
raan = 0;                                    % Right ascension of ascending node, in degrees
argOfPeriapsis = 0;                          % In degrees
trueAnomaly = 343.9391;                      % In degrees
meoSat = satellite(sc, semiMajorAxis, ...
    eccentricity, ...
    inclination, ...
    raan, ...
    argOfPeriapsis, ...
    trueAnomaly, ...
    Name = "MEO Satellite", ...
    OrbitPropagator = "two-body-keplerian");

% Add Interfering Satellite Constellation
interferingSat = satellite(sc,"leoSatelliteConstellation.tle");

% Add Transmitter to MEO Satellites
txMEOFreq = 3e9;                   % In Hz
txMEOSat = transmitter(meoSat, ...
    Frequency = txMEOFreq, ...     % In Hz
    Power = 11);                   % In dBW
gaussianAntenna(txMEOSat, ...
    DishDiameter = 1);             % In m

% Add Transmitter to LEO Satellites
interferenceFreq = 2.99e9;                              % In Hz
rng("default");
txInterferingSat = transmitter(interferingSat, ...
    Frequency = interferenceFreq, ...                   % In Hz
    Power = 10+10*rand(1,numel(interferingSat)));       % In dBW
gaussianAntenna(txInterferingSat, ...
    DishDiameter = 0.2);    

% Add Ground Station
% Add GS(TW_TASA 國家太空中心)
gs = groundStation(sc, ...
    24.4703, ...                    % Latitude in degrees
    121.0015, ...                  % Longitude in degrees
    Name = "TW_GS");

groundStationAntennaType = "Gaussian Antenna";

% Add Receiver to Ground Station
switch groundStationAntennaType
    case {"Gaussian Antenna","Parabolic Reflector"}
        % When Gaussian Antenna or Parabolic Reflector is selected, attach
        % a gimbal to the ground station.
        gim = gimbal(gs, ...
            MountingLocation = [0;0;-5], ... % In m
            MountingAngles = [0;180;0]);     % In degrees

        % Set the gimbal to track the MEO satellite.
        pointAt(gim,meoSat);

        if groundStationAntennaType == "Gaussian Antenna"
            % When Gaussian Antenna is selected

            % Create the receiver object and add it to the gimbal
            rxGs = receiver(gim, ...
                MountingLocation = [0;0;1]); % In m
    
            % Provide the Gaussian Antenna specifications
            gaussianAntenna(rxGs, ...
                DishDiameter = 0.8);         % In m
        else
            % When Parabolic Reflector is selected

            % Size the antenna based on the frequency
            % Requires Antenna Toolbox(TM)
            ant = design(reflectorParabolic,txMEOFreq);
    
            % Create the receiver object and add it to the gimbal
            rxGs = receiver(gim, ...
                Antenna = ant, ...
                MountingLocation = [0;0;1]); % In m
        end
    case "Uniform Rectangular Array"
        % When Uniform Rectangular Array is selected

        % Determine the wavelength of the downlink signal
        c = physconst('LightSpeed');
        lambda = c/txMEOFreq;
         
        % Define array size
        nrow = 8;
        ncol = 8;
         
        % Define element spacing
        drow = lambda/2;
        dcol = lambda/2;
         
        % Create a back-baffled 6-by-6 antenna array
        % Requires Phased Array System Toolbox(TM)
        ant = phased.URA(Size = [nrow ncol], ...
            ElementSpacing = [drow dcol]);
        ant.Element.BackBaffled = true;
        
        % Create the receiver object and add it to the ground station
        rxGs = receiver(gs, ...
            Antenna = ant, ...
            MountingAngles = [0;90;0]); % In degrees
end

% Create Access Analysis Between Interfering Satellite Constellation and Ground Station
ac = access(interferingSat,gs);
ac.LineColor = [1 1 0];              % Yellow

% Set Tracking Targets for Satellites
pointAt([meoSat interferingSat],gs);

% Calculate Weights of Uniform Rectangular Array
if groundStationAntennaType == "Uniform Rectangular Array"
    % Find the LEO satellites that are in the line of sight of the ground
    % station. These satellites are the potential interferers.
    currentInterferingSat = interferingSat(accessStatus(ac,sc.StartTime) == true);
    
    % Calculate the direction of the leo satellite with respect to the
    % array. This is the lookout direction.
    [azd,eld] = aer(rxGs,leoSat,sc.StartTime,CoordinateFrame='body');

    % Calculate the directions of the potentially interfering satellites
    % with respect to the array. These are the null directions.
    [azn,eln] = aer(rxGs,currentInterferingSat,sc.StartTime,CoordinateFrame='body');

    % Calculate the steering vectors for the lookout direction.
    % Requires Phased Array System Toolbox.
    wd = steervec(getElementPosition(ant)/lambda,[wrapTo180(azd);-eld]);

    % Calculate the steering vector for null directions.
    % Requires Phased Array System Toolbox.
    wn = steervec(getElementPosition(ant)/lambda,[wrapTo180(azn)';-eln']);

    % Compute the response of the desired steering at null directions.
    rn = (wn'*wn)\(wn'*wd);

    % Sidelobe canceler - remove the response at null directions.
    w = wd-wn*rn;

    % Assign the weights to the phased array.
    pointAt(rxGs,Weights=w);
end

% Create Desired Downlink
downlink = link(txMEOSat,rxGs);

% Create Interfering Links
lnkInterference = link(txInterferingSat,rxGs);

% Launch Satellite Scenario Viewer
v = satelliteScenarioViewer(sc,ShowDetails=false);

% Visualize Radiation Pattern of Antennas Involved in Downlink
pattern(txMEOSat, ...
    Size = 1000000);        % In m
pattern(rxGs,txMEOFreq, ...
    Size = 1000000);        % In m

% Set Camera to View Ground Station Antenna Radiation Pattern
% Set camera position and orientation to view the ground station antenna
% radiation pattern.
campos(v,-8,172,2500000);
camheading(v,40);
campitch(v,-60);

% Simulate Scenario and Visualize
if groundStationAntennaType == "Gaussian Antenna" || groundStationAntennaType == "Parabolic Reflector"
    play(sc);
    campos(v,-8,172,2500000);
    camheading(v,40);
    campitch(v,-60);
end

if groundStationAntennaType == "Uniform Rectangular Array"
    % Set AutoSimulate to false.
    sc.AutoSimulate = false;

    % Manually step through the simulation.
    while advance(sc)
        % Determine the access status history for each LEO satellite
        % corresponding to the current SimulationTime.
        acStatusHistory = accessStatus(ac);
        acStatus = acStatusHistory(:,end);

        % Determine the LEO satellites that are visible to the ground
        % station. These are the satellites that will potentially
        % interfere with the ground station at the current simulation
        % time.
        currentInterferingSat = interferingSat(acStatus == true);

        % Determine the direction of the leo satellite in the body frame of
        % the Uniform Rectangular Array. This is the lookout direction of
        % the array.
        [azdHistory,eldHistory] = aer(rxGs,leoSat,CoordinateFrame='body');
        azd = azdHistory(:,end);
        eld = eldHistory(:,end);

        % Determine the direction of these interfering satellites in
        % the body frame of the Uniform Rectangular Array. These are
        % the directions in which the array must point a null.
        [aznHistory,elnHistory] = aer(rxGs,currentInterferingSat,CoordinateFrame='body');
        azn = aznHistory(:,end);
        eln = elnHistory(:,end);

        % Calculate the steering vectors for lookout direction.
        % Requires Phased Array System Toolbox.
        wd = steervec(getElementPosition(ant)/lambda,[wrapTo180(azd);-eld]);

        % Calculate the steering vector for null directions.
        % Requires Phased Array System Toolbox.
        wn = steervec(getElementPosition(ant)/lambda,[wrapTo180(azn)';-eln']);

        % Compute the response of desired steering at null direction.
        rn = (wn'*wn)\(wn'*wd);

        % Sidelobe canceler - remove the response at null direction.
        w = wd-wn*rn;

        % Assign the weights to the phased array.
        pointAt(rxGs,Weights=w);
    end
end

% Plot Downlink Closure Status Neglecting Interference
[downlinkStatus,t] = linkStatus(downlink);
plot(t,downlinkStatus,"-g",LineWidth=2);
xlabel("Time");
ylabel("Downlink Closure Status");
title("Link Status as a Function of Time");
grid on;
saveas(gcf,'antenna1.png')

% Calculate Downlink Closure Status with Interference
% Calculate the power at receiver input corresponding to the downlink from
% the leo satellite.
[~,downlinkPowerRxInput] = sigstrength(downlink); % In dBW

% Calculate the interference power at receiver input corresponding to each
% LEO satellite.
[~,interferencePowerRxInput] = sigstrength(lnkInterference); % In dBW

interferencePowerRxInputW = 10.^(interferencePowerRxInput/10); % W
interferencePowerRxInputSumW = sum(interferencePowerRxInputW); % W

txBandwidth = 30e6;                                               % In Hz
interferenceBandWidth = 20e6;                                     % In Hz

% Get the overlap portion of the interfering bandwidth and the bandwidth of
% interest. The assumption is to the have same interference power across
% the whole bandwidth.
overlapFactor = getOverlapFactor(txMEOFreq,txBandwidth, ...
    interferenceFreq,interferenceBandWidth);

% Get the interference power that contributes as interference to the signal
% of interest from the total interference power
interferencePowerRxInputActual = interferencePowerRxInputSumW*overlapFactor; % In W

% Calculate the thermal noise at the ground station receiver input.
T = HelperGetNoiseTemperature(txMEOFreq,rxGs); % In K
kb = physconst("Boltzmann");
thermalNoise = kb*T*txBandwidth;               % In W

% Calculate the noise plus interference power at the receiver input.
noisePlusInterferencePowerW = thermalNoise + interferencePowerRxInputActual; % In W
noisePlusInterferencePower = 10*log10(noisePlusInterferencePowerW);          % In dBW

% Calculate loss that occurs between the receiver input and the demodulator
% input.
rxGsLoss = rxGs.SystemLoss - rxGs.PreReceiverLoss;

% Calculate C/(N0+I0) at the demodulator input.
CNoPlusInterference = downlinkPowerRxInput - ...
    noisePlusInterferencePower + 10*log10(txBandwidth) - rxGsLoss;

bitRate = txMEOSat.BitRate;
ebNoPlusInterference = CNoPlusInterference - 10*log10(bitRate) - 60;

marginWithInterference = ebNoPlusInterference - rxGs.RequiredEbNo;

downlinkStatusWithInterference = marginWithInterference >= 0;

% Calculate Energy per Bit to Noise Power Spectral Density Ratio
ebnoDownlink = ebno(downlink);            % In dB
ebnoInterference = ebno(lnkInterference); % In dB

% Plot Downlink Closure Status with Interference
plot(t,downlinkStatusWithInterference,"-r",t,downlinkStatus,"--g",LineWidth=2);
legend("Interference accounted","Interference neglected");
xlabel("Time");
ylabel("Downlink Closure Status");
title("Link Status as a Function of Time");
ylim([0 1.2]);
grid on
saveas(gcf,'antenna2.png')

if groundStationAntennaType == "Gaussian Antenna" || groundStationAntennaType == "Parabolic Reflector"
    v.CurrentTime = datetime(2024,4,18,22,0,0);
end

if groundStationAntennaType == "Uniform Rectangular Array"
    % Set AutoSimulate to true.
    sc.AutoSimulate = true;

    % Set viewer CurrentTime to 10:00 PM.
    time = datetime(2024,4,18,22,0,0);
    v.CurrentTime = time;

    % Calculate the weights and assign them to the array.
    currentInterferingSat = interferingSat(accessStatus(ac,time) == true);
    [azd,eld] = aer(rxGs,meoSat,time,CoordinateFrame='body');
    [azn,eln] = aer(rxGs,currentInterferingSat,time,CoordinateFrame='body');

    % Requires Phased Array System Toolbox.
    wd = steervec(getElementPosition(ant)/lambda,[wrapTo180(azd);-eld]);
    wn = steervec(getElementPosition(ant)/lambda,[wrapTo180(azn)';-eln']);
    
    rn = (wn'*wn)\(wn'*wd);
    w = wd-wn*rn;
    pointAt(rxGs,Weights=w);

    % Make the radiation pattern opaque.
    pattern(rxGs,txMEOFreq, ...
        Size = 1000000, ...
        Transparency = 1);
end

% Calculate and Plot Carrier to Noise Ratio and Carrier to Noise Plus Interference Ratio
% Calculate the carrier to noise power spectral density ratio.
CNoDownlink = ebnoDownlink + 10*log10(bitRate) + 60;

% Calculate the carrier to noise ratio.
cByN = CNoDownlink - 10*log10(txBandwidth);
cByNPlusI = CNoPlusInterference - 10*log10(txBandwidth);

plot(t,cByNPlusI,"-r",t,cByN,"--g",LineWidth=2);
legend("CNIR", "CNR",Location="south");
xlabel("Time");
ylabel("CNR or CNIR (dB)");
title("CNR and CNIR vs. Time for " + groundStationAntennaType);
grid on
saveas(gcf,'antenna3.png')

maxInterferencePowerRxInput = max(interferencePowerRxInput,[],'all');
disp("The maximum power at ground station receiver input from the interfering LEO satellites over the entire scenario duration is " + maxInterferencePowerRxInput + " dBW.");

% Compare Link Margins with and without Interference
marginWithoutInterference = ebnoDownlink - rxGs.RequiredEbNo;

plot(t,marginWithInterference,"-r",t,marginWithoutInterference,"--g",LineWidth=2);
legend("With interference","Without interference",Location="south");
xlabel("Time");
ylabel("Margin (dB)");
title("Link Margin vs. Time for " + groundStationAntennaType);
grid on
saveas(gcf,'antenna4.png')

% Helper Functions
function overlapFactor = getOverlapFactor(txFreq,txBW,interferenceFreq,interferenceBW)
% getOverlapFactor provides the amount of interference bandwidth overlapped
% with transmission bandwidth

    txFreq_Limits = [txFreq-(txBW/2) txFreq+(txBW/2)];
    interferenceFreq_Limits = [interferenceFreq-(interferenceBW/2) ...
        interferenceFreq+(interferenceBW/2)];
    if (interferenceFreq_Limits(2) < txFreq_Limits(1)) || ...
            (interferenceFreq_Limits(1) > txFreq_Limits(2))
        % If no overlap exists between transmission bandwidth and
        % interfering bandwidth, then overlap factor is 0
        overlapFactor = 0;
    elseif (interferenceFreq_Limits(2) <= txFreq_Limits(2)) && ...
            (interferenceFreq_Limits(1) >= txFreq_Limits(1))
        % If interfering bandwidth lies completely within transmission
        % bandwidth, then overlap factor is 1
        overlapFactor = 1;
    elseif (interferenceFreq_Limits(2) > txFreq_Limits(2)) && ...
            (interferenceFreq_Limits(1) < txFreq_Limits(1))
        % If transmission bandwidth lies completely within interfering
        % bandwidth, then overlap factor is the ratio of transmission
        % bandwidth with that of interference bandwidth
        overlapFactor = txBW/interferenceBW;
    elseif (interferenceFreq_Limits(2) <= txFreq_Limits(2)) && ...
            (interferenceFreq_Limits(1) <= txFreq_Limits(1))
        % If the start edge of transmission bandwidth lies within
        % interfering bandwidth, then overlap factor is the ratio of
        % difference from last edge of interfering bandwidth and first edge
        % of signal bandwidth, with that of interference bandwidth
        overlapFactor = (interferenceFreq_Limits(2)-txFreq_Limits(1))/interferenceBW;
    else
        % If the last edge of transmission bandwidth lies within
        % interfering bandwidth, then overlap factor is the ratio of difference
        % from last edge of signal bandwidth and first edge of interfering
        % bandwidth, with that of interference bandwidth
        overlapFactor = (-interferenceFreq_Limits(1)+txFreq_Limits(2))/interferenceBW;
    end

end