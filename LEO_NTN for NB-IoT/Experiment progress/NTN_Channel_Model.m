% Set NTN Channel Common Parameters
commonParams = struct;
commonParams.CarrierFrequency = 2e9;              % In Hz
commonParams.ElevationAngle = 50;                 % In degrees
commonParams.SatelliteAltitude = 600000;          % In m
commonParams.SatelliteSpeed = 7562.2;             % In m/s
commonParams.MobileSpeed = 3*1000/3600;           % In m/s
commonParams.SampleRate = 7680000;                % In Hz
% Set the random stream and seed, for reproducibility
commonParams.RandomStream = "mt19937ar with seed";
commonParams.Seed = 73;
% Set the number of sinusoids used in generation of Doppler spread
commonParams.NumSinusoids = 48;

% Set NTN Narrowband Channel Parameters
% Initialize the NTN flat fading narrowband channel parameters in a
% structure
ntnNarrowbandParams = commonParams;
ntnNarrowbandParams.NTNChannelType = "Narrowband";
ntnNarrowbandParams.Environment = "Urban";
ntnNarrowbandParams.AzimuthOrientation = 0;
ntnNarrowbandParams.FadingTechnique = "Sum of sinusoids";

% Set the below parameters when Environment is set to Custom
ntnNarrowbandParams.StateDistribution = [3.0639 2.9108; 1.6980 1.2602];
ntnNarrowbandParams.MinStateDuration = [10 6];
ntnNarrowbandParams.DirectPathDistribution = [-1.8225 -15.4844; 1.1317 3.3245];
ntnNarrowbandParams.MultipathPowerCoefficients = [-0.0481 0.9434; -14.7450 -1.7555];
ntnNarrowbandParams.StandardDeviationCoefficients = [-0.4643 -0.0798; 0.3334 2.8101];
ntnNarrowbandParams.DirectPathCorrelationDistance = [1.7910 1.7910];
ntnNarrowbandParams.TransitionLengthCoefficients = [0.0744; 2.1423];
ntnNarrowbandParams.StateProbabilityRange = [0.05 0.1; 0.95 0.9];

ntnNarrowbandChan = HelperSetupNTNChannel(ntnNarrowbandParams);

p681ChannelInfo = info(ntnNarrowbandChan.BaseChannel);

% Generate NTN Narrowband Channel
% Generate a random input
rng(commonParams.Seed);
in = complex(randn(commonParams.SampleRate,1), ...
    randn(commonParams.SampleRate,1));

[narrowbandOut,narrowbandPathGains,narrowbandSampleTimes] = ...
    HelperGenerateNTNChannel(ntnNarrowbandChan,in);

% Visualize NTN Narrowband Channel Received Spectrum
ntnNarrowbandAnalyzer = spectrumAnalyzer( ...
    SampleRate = ntnNarrowbandParams.SampleRate);
ntnNarrowbandAnalyzer.Title = "Received Signal Spectrum " ...
    + ntnNarrowbandChan.ChannelName;
ntnNarrowbandAnalyzer.ShowLegend = true;
ntnNarrowbandAnalyzer.ChannelNames = "Rx Antenna 1";
ntnNarrowbandAnalyzer(narrowbandOut)

%% NTN TDL Channel
% This example supports all the four channel profiles of the NTN TDL channel, defined in 3GPP TR 38.811 Section 6.9.2
% These four channel profiles are: NTN-TDL-A ~ D

%% Set NTN TDL Channel Parameters
% Initialize the NTN TDL channel parameters in a structure
ntnTDLParams = commonParams;
ntnTDLParams.NTNChannelType = "TDL";
ntnTDLParams.DelayProfile = "NTN-TDL-D";
ntnTDLParams.DelaySpread = 30e-9;                          % In s
ntnTDLParams.TransmissionDirection = "Downlink";
ntnTDLParams.MIMOCorrelation = "Low";
ntnTDLParams.Polarization = "Co-Polar";
% Modify the below parameters, when DelayProfile is set to Custom
ntnTDLParams.PathDelays = 0;                               % In s
ntnTDLParams.AveragePathGains = 0;                         % In dB
ntnTDLParams.FadingDistribution = "Rayleigh";

% Set the antenna configuration
% Modify the below parameters, when MIMOCorrelation is set to a value other
% than Custom
ntnTDLParams.NumTransmitAntennas = 1;
ntnTDLParams.NumReceiveAntennas = 2;
% Modify the below parameters, when MIMOCorrelation is set to Custom and
% Polarization is set to Co-Polar or Cross-Polar
ntnTDLParams.TransmitCorrelationMatrix = 1;
ntnTDLParams.ReceiveCorrelationMatrix = [1 0; 0 1];
% Modify the below parameters, when MIMOCorrelation is set to Custom and
% Polarization is set to Cross-Polar
ntnTDLParams.TransmitPolarizationAngles = [45 -45];        % In degrees
ntnTDLParams.ReceivePolarizationAngles = [90 0];           % In degrees
ntnTDLParams.XPR = 10;                                     % In dB
% Modify the below parameters, when both MIMOCorrelation and Polarization
% are set to Custom
ntnTDLParams.SpatialCorrelationMatrix = [1 0; 0 1];

ntnTDLChan = HelperSetupNTNChannel(ntnTDLParams);

tdlChanInfo = info(ntnTDLChan.BaseChannel);

%% Generate NTN TDL Channel
% Generate a random input
rng(commonParams.Seed);
in = complex(randn(commonParams.SampleRate,tdlChanInfo.NumTransmitAntennas), ...
    randn(commonParams.SampleRate,tdlChanInfo.NumTransmitAntennas));
% Generate the faded waveform for NTN TDL channel
[tdlOut,tdlPathGains,tdlSampleTimes] = HelperGenerateNTNChannel(ntnTDLChan,in);

%% Visualize NTN TDL Channel Received Spectrum
ntnTDLAnalyzer = spectrumAnalyzer(SampleRate = ntnTDLParams.SampleRate);
ntnTDLAnalyzer.Title = "Received Signal Spectrum " ...
    + ntnTDLChan.ChannelName;
ntnTDLAnalyzer.ShowLegend = true;
for nRx = 1:size(tdlOut,2)
    ntnTDLAnalyzer.ChannelNames{nRx} = "Rx Antenna " + nRx;
end
ntnTDLAnalyzer(tdlOut)






