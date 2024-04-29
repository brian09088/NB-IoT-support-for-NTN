clc,clear,close all
%% create satellite scenario
startTime = datetime(2024,3,16,8,0,0);              
stopTime = startTime + minutes(10);                  
sampleTime = 60;   % seconds                 
sc = satelliteScenario(startTime,stopTime,sampleTime);

%% Add Medium-Earth Orbit Satellite
Re = 6378e3; % Radius of Earth
LEO_altitude = 650e3;
semiMajorAxis = Re + LEO_altitude;    % In m
eccentricity = 0.000258;
inclination = 51.65;          % In degrees
raan = 0;                  % Right ascension of ascending node, in degrees
argOfPeriapsis = 0;        % In degrees
trueAnomaly = 0;           % In degrees
leoSat = satellite(sc, semiMajorAxis, ...
    eccentricity, ...
    inclination, ...
    raan, ...
    argOfPeriapsis, ...
    trueAnomaly, ...
    Name = "LEO Satellite", ...
    OrbitPropagator = "two-body-keplerian");

%% Add Interfering Satellite Constellation
interferingSat = satellite(sc,"leoSatelliteConstellation.tle");

%% Add Transmitter to LEO Satellite
interferenceFreq = 2.99e9;                              % In Hz
rng("default");
txInterferingSat = transmitter(interferingSat, ...
    Frequency = interferenceFreq, ...                   % In Hz
    Power = 10+10*rand(1,numel(interferingSat)));       % In dBW
gaussianAntenna(txInterferingSat, ...
    DishDiameter = 0.2);                                % In m

%% Add GS(TW_TASA 國家太空中心)
gs = groundStation(sc, ...
    24.4703, ...                    % Latitude in degrees
    121.0015, ...                  % Longitude in degrees
    Name = "TW_GS");

%% groundStationAntennaType
groundStationAntennaType = "Uniform Rectangular Array";

%% Add receiver to GS
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
            ant = design(reflectorParabolic,interferenceFreq);
    
            % Create the receiver object and add it to the gimbal
            rxGs = receiver(gim, ...
                        Antenna = ant, ...
                        MountingLocation = [0;0;1]); % In m
        end
    case "Uniform Rectangular Array"
        % When Uniform Rectangular Array is selected

        % Determine the wavelength of the downlink signal
        c = physconst('LightSpeed');
        lambda = c/interferenceFreq;
         
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

c = physconst('LightSpeed');
lambda = c/interferenceFreq;
 
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

% Create Access Analysis Between Interfering Satellite Constellation and Ground Station
ac = access(interferingSat,gs);
ac.LineColor = [1 1 0];

% set tracking targets for satellite
pointAt([leoSat interferingSat],gs);

downlink = link(txInterferingSat,rxGs);

lnkInterference = link(txInterferingSat,rxGs);

v = satelliteScenarioViewer(sc,ShowDetails=false);

pattern(txInterferingSat, ...
    Size = 1000000);        % In m
pattern(rxGs,interferenceFreq, ...
    Size = 1000000);        % In m

% Set camera position and orientation to view the ground station antenna
% radiation pattern.
campos(v,-8,172,2500000);
camheading(v,40);
campitch(v,-60);

if groundStationAntennaType == "Gaussian Antenna" || groundStationAntennaType == "Parabolic Reflector"
    play(sc);
    campos(v,-8,172,2500000);
    camheading(v,40);
    campitch(v,-60);
end

if groundStationAntennaType == "Gaussian Antenna" || groundStationAntennaType == "Parabolic Reflector"
    v.CurrentTime = datetime(2024,4,18,22,0,0);
end

if groundStationAntennaType == "Uniform Rectangular Array"
    play(sc);
    campos(v,-8,172,2500000);
    camheading(v,40);
    campitch(v,-60);
end

%% play scenario video
% can be modify to other type
if groundStationAntennaType == "Uniform Rectangular Array"
    % Set AutoSimulate to false.
    sc.AutoSimulate = false;

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
    
        % Determine the direction of the MEO satellite in the body frame of
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