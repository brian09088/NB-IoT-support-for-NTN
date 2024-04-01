clc,clear,close all

satelliteParamsSource = "Set 3 LEO-600";
%% Satellite parameters
% When you set satelliteParamsSource to Custom, provide these fields.
satellite = struct;
satellite.EIRPDensity = 59;   % dBW/MHz
satellite.RxGByT = 19;        % dB/K
satellite.Altitude = 600e3; % 60000m = 600km

%% Set User Equipment Payload Characteristics
% By default, this example uses the UE power class PC3 with a noise figure of 7 dB
ue = struct;
ue.TxPower = 23;               % dBm
ue.TxGain = 0;                 % dBi
ue.TxCableLoss = 0;            % dB
ue.RxNoiseFigure = 7;          % dB
ue.RxGain = 0;                 % dBi
ue.RxAntennaTemperature = 290; % K
ue.RxAmbientTemperature = 290; % K

% Set Link Characteristics
% Set the link characteristics by defining the link direction, elevation angle, bandwidth, frequency, and atmospheric losses.
% The transmission bandwidths for NB-IoT are:
% Downlink — 180 kHz
% Uplink — 3.75 kHz, 15 kHz, 45 kHz, 90 kHz, or 180 kHz

% Enable or disable the ITU-R P.618 propagation losses

useP618PropagationLosses = false;           % 0 (false), 1 (true)

%% Set the link characteristics
link = struct;
link.Direction = "downlink";    % "uplink", "downlink"
link.ElevationAngle = [10.95 20];         % degrees (Vector of elevation angles)
link.Frequency = 2e9;                     % Hz
link.Bandwidth = 180e3;                   % Hz
link.ShadowMargin = 3;                    % dB
link.AdditionalLosses = 0;                % dB
link.PolarizationLoss = 3;                % dB

if useP618PropagationLosses == 0
    % Set these fields when you set useP618PropagationLosses to false.
    link.ScintillationLosses = 2.2;           % dB
    link.AtmosphericLosses = 0.2;             % dB
else
    % Set these fields when you set useP618PropagationLosses to true. In
    % this case, the example uses digital maps to calculate the
    % attenuation.
    link.P618Configuration = p618Config;
    link.P618Configuration.Latitude = 51.5;                   % degrees
    link.P618Configuration.Longitude = -0.14;                 % degrees
    link.P618Configuration.GasAnnualExceedance = 1;
    link.P618Configuration.CloudAnnualExceedance = 1;
    link.P618Configuration.ScintillationAnnualExceedance = 1;
    link.P618Configuration.TotalAnnualExceedance = 1;
    link.P618Configuration.PolarizationTiltAngle = 0;         % degrees
    link.P618Configuration.AntennaDiameter = 1;               % meters
    link.P618Configuration.AntennaEfficiency = 0.5;
end

%% Set NB-IoT Data Channel Parameters
% To calculate the reference carrier-to-noise ratio (CNR) for an NB-IoT data channel transmission
% set these NB-IoT parameters
% - Modulation scheme
% - Transport block size (Nbits)
% - Number of symbols used for transmission (Nsymbols)
% - Number of repetitions (NRep)
% - Number of subframes in downlink (NSF)
% - Number of resource units in uplink (NRU)

modulation = "QPSK"; % Modulation scheme
nBits = 208;                       % Number of useful bits (transport block size)
nSymbols = 160;                    % Number of symbols used for transmission
nRep = 1;                          % Number of repetitions
nSF = 8;                           % Number of subframes (used in downlink)
nRU = 1;                           % Number of resource units (used in uplink)

%% Set these NB-IoT waveform parameters.
% - Number of used data subcarriers (NDSC)
% - Number of fast Fourier transform (FFT) bins (NFFT)
% - Data symbol duration (Td)
% - Cyclic prefix duration (TCP)
% - Oversampling factor (OSR)

nDSC = 72;                         % Number of data subcarriers (6 resource blocks)
nFFT = 128;                        % FFT length
tsamp = 1/1.92e6;                  % Sample time (in s)
td = nFFT*tsamp;                   % Data symbol duration (in s)
nCP = 9;                           % Number of cyclic prefix samples
tCP = nCP*tsamp;                   % Cyclic prefix duration (in s)
osr = 1;                           % Oversampling factor

% In the presence of additive white Gaussian noise without any channel coding
% the modulation schemes return these results:
% - BPSK achieves a bit error rate of 1e-6 at an Eb/No of 10.5dB
% - QPSK achieves a bit error rate of 1e-6 at an Eb/No of 10.5dB
% - 16-QAM achieves a bit error rate of 1e-6 at an Eb/No of 14.4dB

%% Reference Eb/No in dB
if modulation == "16-QAM"
    % Modulation scheme is 16-QAM
    ebnoRef = 14.4;
else 
    % Modulation scheme is either BPSK or QPSK
    ebnoRef = 10.5;
end

%% Compute Reference CNR
% To compute the reference CNR in dB, use this equation:
% (C/N)ref = (Eb/No)ref + 10*log(m*Reff) + 10*log(Ndsc/Nfft) + 10*log(Td/(Td+Tcp)) - 10*log(OSR)

% (Eb/No)ref is the value in dB for the desired bit error rate
%　m is the modulation order
% Reff is the effective code rate

%% The effective code rate for NB-IoT downlink:
% Reff = (nBits + nCRC)/(nSF*nSymbols*m*nRep)
% - Nbits is the number of useful bits to transmit (transport block size)
% - NCRC is the number of bits used for check code (24 for NB-IoT)
% - NSF is the number of allocated subframes
% - Nsymbols is the number of symbols
% - Nrep is the number of repetitions

%% The effective code rate for NB-IoT uplink:
% Reff = (nBits + nCRC)/(nRU*nSymbols*m*nRep)
% NRU is the number of allocated resource units
% For 1-tone transmission, Nsymbols is 96
% for 3-tone transmission, Nsymbols is 144

% Modulation order
m = 1;
if modulation == "QPSK"
    m = 2;
elseif modulation == "16-QAM"
    m = 4;
end

% Number of CRC bits
nCRC = 24;

% Calculate the effective code rate and the reference CNR
if link.Direction == "downlink"
    Reff = (nBits + nCRC)/(nSF*nSymbols*m*nRep);
else % "uplink"
    Reff = (nBits + nCRC)/(nRU*nSymbols*m*nRep);
end
cnrRef = ebnoRef + 10*log10(m*Reff) + 10*log10(nDSC/nFFT) + ...
        10*log10(td/(td + tCP)) - 10*log10(osr);

disp("Reference CNR: " + cnrRef + " dB")

%% Compute Link Budget
% EIRP is the effective isotropic radiated power in dBW.
% G/T is the antenna-gain-to-noise temperature in dB/K.
% k is the Boltzmann constant with the value of -228.6 dBW/K/Hz.
% PL_FS is the free space path loss (FSPL) in dB.
% PL_A is the atmospheric path loss due to gases and shadowing margin in dB.
% PL_SM is the shadowing margin in dB.
% PL_SL is the scintillation loss in dB.
% PL_AD is the additional loss in dB.
% B is the channel bandwidth in dBHz.

% To calculate antenna-gain-to-noise-temperature
% GR is the receive antenna gain in dBi.
% Nf is the noise figure in dB.
% T0 is the ambient temperature in degrees Kelvin.
% Ta is the antenna temperature in degrees Kelvin.

% The receive antenna gain depends on the type of antenna used, and accounts for polarization loss.
% EIRP = PT - LC + GT
% PT is the transmit antenna power in dBW.
% LC is the cable loss in dB.
% GT is the transmit antenna gain in dBi.

% Get the satellite parameters based on satelliteParamsSource
if satelliteParamsSource ~= "Custom"
    satellite = getSatelliteParams(satelliteParamsSource);
end

% Calculate EIRP when the satellite is a transmitter (dB)
% To get a value in decibels from density: 
% Value in dB = Value in dB/MHz + 10*log10[BW in MHz]
satellite.EIRP = satellite.EIRPDensity + 10*log10(link.Bandwidth/1e6);

% Calculate EIRP when the UE is a transmitter (dB)
% To get a value in decibels from dBm: dB = dBm - 30
ue.EIRP = (ue.TxPower-30) + ue.TxGain - ue.TxCableLoss;

% Calculate the gain-to-noise temperature or the figure of merit for the UE
% as a receiver.
ue.RxGByT = ue.RxGain - ue.RxNoiseFigure ...
    - 10*log10(ue.RxAmbientTemperature + ...
    (ue.RxAntennaTemperature-ue.RxAmbientTemperature)*10^(-0.1*ue.RxNoiseFigure));

% Set the transmitter and the receiver based on the link direction
if link.Direction == "uplink"
    tx = ue;
    rx = satellite;
else
    tx = satellite;
    rx = ue;
end

% Range of elevation angles
elevAngles = link.ElevationAngle(:);
numElevAngles = numel(elevAngles);

% Calculate the distance from the satellite to the ground station for all
% elevation angles
re = physconst("earthradius");
c = physconst("lightspeed");
h = satellite.Altitude;
d = -re*sind(elevAngles) + sqrt((re*sind(elevAngles)).^2 + h*h + 2*re*h);

% Get the total atmospheric losses for each elevation angle
totalAtmosphericLoss = zeros(numElevAngles,1);
if useP618PropagationLosses == 1
    % Download and extract the digital maps, if not available on the path
    maps = exist("maps.mat","file");
    p836 = exist("p836.mat","file");
    p837 = exist("p837.mat","file");
    p840 = exist("p840.mat","file");
    matFiles = [maps p836 p837 p840];
    if ~all(matFiles)
        if ~exist("ITURDigitalMaps.tar.gz","file")
            url = "https://www.mathworks.com/supportfiles/spc/P618/ITURDigitalMaps.tar.gz";
            websave("ITURDigitalMaps.tar.gz",url);
            untar("ITURDigitalMaps.tar.gz")
        else
            untar("ITURDigitalMaps.tar.gz")
        end
    end
    link.P618Configuration.Frequency = link.Frequency;
    elevAnglesToConsider = elevAngles;
    if any(elevAngles < 5)
        warning("The prediction method for scintillation losses is valid for elevation " + ... 
            "angle greater than 5 degree. For elevation angle less than 5 degree, the " + ... 
            "nearest valid value of 5 degree will be used in the computation.")
        elevAnglesToConsider(elevAngles < 5) = 5;
    end
    for index = 1:numElevAngles
        link.P618Configuration.ElevationAngle = elevAnglesToConsider(index);
        pl = p618PropagationLosses(link.P618Configuration);
        totalAtmosphericLoss(index) = pl.At;
    end
else
    totalAtmosphericLoss(:) = link.AtmosphericLosses + link.ScintillationLosses;
end

% Define the configuration parameters
config = satelliteCNRConfig;
config.TransmitterPower = tx.EIRP;
config.TransmitterAntennaGain = 0;
config.Frequency = link.Frequency/1e9;          % GHz
config.GainToNoiseTemperatureRatio = rx.RxGByT;
config.Bandwidth = link.Bandwidth/1e6;          % MHz

% Compute the CNR and the free space path loss for each elevation angle
cnr = zeros(numElevAngles,1);
pathLoss = cnr;
for index = 1:numElevAngles
    config.Distance = d(index)/1e3;                                         % km
    config.MiscellaneousLoss = totalAtmosphericLoss(index) + ...
         link.PolarizationLoss + link.ShadowMargin + link.AdditionalLosses;
    [cnr(index),cnrInfo] = satelliteCNR(config);
    pathLoss(index) = cnrInfo.FSPL;
end

% Report the results in a table
table(elevAngles,cnr,pathLoss,VariableNames=["Elevation Angle (degrees)", "CNR (dB)", "FSPL (dB)"])

% Plot CNR as a function of the elevation angle
plot(elevAngles,cnr,"*")
title("CNR As a Function of Elevation Angle")
xlabel("Elevation Angle (degrees)")
ylabel("CNR (dB)")
ylim([min(cnr)-0.2 max(cnr)+0.2])
xlim([min(elevAngles)-1 max(elevAngles)+1])
grid on

%% Compute Link Margin and NB-IoT Repetitions
% The link margin (in dB) is the difference between the calculated CNR and the reference CNR. 
% The link is closed when the link margin is zero or positive.
% To calculate the link margin, use this equation: 
% LinkMargin = C/N-(C/N)ref

% Compute link margin
linkMargin = cnr - cnrRef;

%% Minimum number of additional repetitions required for link closure
minRepetitions = 10.^(-linkMargin./10);
idx = linkMargin >= 0;
additionalRepetitions = minRepetitions;
additionalRepetitions(idx) = 0;
% When the link margin is negative, improve the CNR by adding repetitions.
% Calculate the required number of additional repetitions.
additionalRepetitions(~idx) = ceil(nRep*(additionalRepetitions(~idx)-1));

% Report the results in a table with this format.
% Elevation Angle | Link Margin | Min. Additional Repetitions Required
table(elevAngles,linkMargin,additionalRepetitions, ...
    VariableNames=["Elevation Angle (degrees)", ...
    "Link Margin (dB)", "NRep_Add"])

%% Local Function
% getSatelliteParams — Gets the satellite parameters defined for the selected satelliteParamsSource
function satellite = getSatelliteParams(satelliteParamsSource)
    
    if satelliteParamsSource == "Set 1 GEO"
        % Table 6.2-4, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 59;    % dBW/MHz
        satellite.RxGByT = 19;         % dB/K
        satellite.Altitude = 35786e3;  % m
    elseif satelliteParamsSource == "Set 1 LEO-1200"
        % Table 6.2-4, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 40;    % dBW/MHz
        satellite.RxGByT = 1.1;        % dB/K
        satellite.Altitude = 1200e3;   % m
    elseif satelliteParamsSource == "Set 1 LEO-600"
        % Table 6.2-4, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 34;    % dBW/MHz
        satellite.RxGByT = 1.1;        % dB/K
        satellite.Altitude = 600e3;    % m
    elseif satelliteParamsSource == "Set 2 GEO"
        % Table 6.2-5, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 53.5;  % dBW/MHz
        satellite.RxGByT = 14;         % dB/K
        satellite.Altitude = 35786e3;  % m
    elseif satelliteParamsSource == "Set 2 LEO-1200"
        % Table 6.2-5, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 34;    % dBW/MHz
        satellite.RxGByT = -4.9;       % dB/K
        satellite.Altitude = 1200e3;   % m
    elseif satelliteParamsSource == "Set 2 LEO-600"
        % Table 6.2-5, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 28;    % dBW/MHz
        satellite.RxGByT = -4.9;       % dB/K
        satellite.Altitude = 600e3;    % m
    elseif satelliteParamsSource == "Set 3 GEO"
        % Table 6.2-6, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 59.8;  % dBW/MHz
        satellite.RxGByT = 16.7;       % dB/K
        satellite.Altitude = 35786e3;  % m
    elseif satelliteParamsSource == "Set 3 LEO-1200"
        % Table 6.2-6, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 33.7;  % dBW/MHz
        satellite.RxGByT = -12.8;      % dB/K
        satellite.Altitude = 1200e3;   % m
    elseif satelliteParamsSource == "Set 3 LEO-600"
        % Table 6.2-6, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 28.3;  % dBW/MHz
        satellite.RxGByT = -12.8;      % dB/K
        satellite.Altitude = 600e3;    % m
    elseif satelliteParamsSource == "Set 4 LEO-600"
        % Table 6.2-7, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 21.45; % dBW/MHz
        satellite.RxGByT = -18.6;      % dB/K
        satellite.Altitude = 600e3;    % m
    else % "Set 5 MEO-10000"
        % Table 6.2-8, TR 36.763
        satellite = struct;
        satellite.EIRPDensity = 45.4;  % dBW/MHz
        satellite.RxGByT = 3.8;        % dB/K
        satellite.Altitude = 10000e3;  % m
    end

end




