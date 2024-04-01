function ntnChan = HelperSetupNTNChannel(chanParams)
%HelperSetupNTNChannel Set up NTN channel
%
%   Note: This is a helper and its API and/or functionality may change
%   in subsequent releases.
%
%   NTNCHAN = HelperSetupNTNChannel(CHANPARAMS) sets up the New Radio
%   (NR) non-terrestrial network (NTN) narrowband or tapped delay line
%   (TDL) channel, as mentioned in TR 38.811, given the set of channel
%   parameters CHANPARAMS. NTNCHAN is a structure consisting of building
%   blocks of NTN channel.
%
%   NTNCHAN is an output structure with fields:
%   ChannelName           - Name of the NTN channel.
%   BaseChannel           - Base channel model to generate the path gains
%                           due to mobile movement. It is either
%                           p681LMSChannel or nrTDLChannel System
%                           object.
%   ChannelFilter         - Filter signal using multipath gains at
%                           specified path delays using this
%                           comm.ChannelFilter System object.
%   SatelliteDopplerShift - Doppler shift due to satellite movement.
%   MobileDopplerSpread   - Maximum Doppler shift due to mobile or UE
%                           movement.
%
%   CHANPARAMS is an input structure with the fields:
%
%   NTNChannelType    - Type of NTN channel ("Narrowband", "TDL")
%   CarrierFrequency  - Carrier frequency of input signal (Hz)
%   ElevationAngle    - Elevation angle (degrees)
%   SatelliteSpeed    - Relative speed of satellite with respect to earth (m/s)
%   SatelliteAltitude - Satellite altitude from the surface of earth (m)
%   MobileSpeed       - Speed of mobile terminal (m/s)
%   SampleRate        - Optional. Input signal sample rate (Hz)
%                       (default 30.72e6 for TDL, 7.68e6 for narrowband)
%   OutputDataType    - Optional. Data type of channel output
%                       ("double" (default), "single")
%
%   Specify the following fields when NTNChannelType is set to
%   "Narrowband". For more information of the optional fields, look at the
%   help of <a href="matlab:help p681LMSChannel">p681LMSChannel</a>.
%   AzimuthOrientation            - Azimuth orientation (degrees)
%   Environment                   - Type of propagation environment ("Urban",
%                                   "Suburban", "RuralWooded", "Village", "Train",
%                                   "Residential", "Highway", "Rural", "Custom")
%   StateDistribution             - Optional. Parameters of state distribution (dB)
%                                   (default [3.0639 2.9108; 1.6980 1.2602])
%   MinStateDuration              - Optional. Minimum duration of each state (m)
%                                   (default [10 6])
%   DirectPathDistribution        - Optional. Parameters of direct path amplitude
%                                   distribution (dB)
%                                   (default [-1.8225 -15.4844; 1.1317 3.3245])
%   MultipathPowerCoefficients    - Optional. Coefficients to calculate multipath
%                                   power (default [-0.0481 0.9434; -14.7450 -1.7555])
%   StandardDeviationCoefficients - Optional. Coefficients to compute standard
%                                   deviation of direct path amplitude
%                                   (default [-0.4643 -0.0798; 0.3334 2.8101])
%   DirectPathCorrelationDistance - Optional. Direct path amplitude correlation
%                                   distance (m) (default [1.7910 1.7910])
%   TransitionLengthCoefficients  - Optional. Coefficients to compute transition
%                                   length (default [0.0744; 2.1423])
%   StateProbabilityRange         - Optional. Minimum and maximum probability of
%                                   each state (default [0.05 0.1; 0.95 0.9])
%   ChannelFiltering              - Optional. Perform filtering of input signal
%                                   (true (default), false)
%   NumSamples                    - Optional. Number of time samples (default 7680)
%   FadingTechnique               - Optional. Channel model fading technique
%                                   ("Filtered Gaussian noise" (default),
%                                   "Sum of sinusoids")
%   NumSinusoids                  - Optional. Number of sinusoids (default 48)
%   RandomStream                  - Optional. Source of random number stream
%                                   ("Global stream" (default), "mt19937ar with seed")
%   Seed                          - Optional. Initial seed of mt19937ar random
%                                   number stream generator (default 73)
%
%   Specify the following fields when NTNChannelType is set to "TDL". For
%   more information of the optional fields, look at the help of <a
%   href="matlab:help nrTDLChannel">nrTDLChannel</a>.
%   DelayProfile               - NTN TDL delay profile ("NTN-TDL-A",
%                                "NTN-TDL-B", "NTN-TDL-C", "NTN-TDL-D", "Custom")
%   PathDelays                 - Optional. Discrete path delay vector (s)
%                                (default 0)
%   AveragePathGains           - Optional. Average path gain vector (dB)
%                                (default 0)
%   FadingDistribution         - Optional. Rayleigh or Rician fading
%                                ("Rayleigh" (default), "Rician")
%   KFactorFirstTap            - Optional. K-factor of first tap (dB)
%                                (default 13.3 dB)
%   DelaySpread                - Optional. Desired delay spread (s)
%                                (default 30 ns)
%   KFactorScaling             - Optional. Enable K-factor scaling (logical)
%                                (default false)
%   KFactor                    - Optional. Desired Rician K-factor (dB)
%                                (default 10.224)
%   MIMOCorrelation            - Optional. Correlation between UE and BS antennas
%                                ("Low" (default), "Medium", "Medium-A",
%                                "UplinkMedium", "High", "Custom")
%   Polarization               - Optional. Antenna polarization arrangement
%                                ("Co-Polar" (default), "Cross-Polar",
%                                "Custom")
%   TransmissionDirection      - Optional. Transmission direction (Uplink/Downlink)
%                                (default "Downlink")
%   NumTransmitAntennas        - Optional. Number of transmit antennas (default 1)
%   NumReceiveAntennas         - Optional. Number of receive antennas (default 2)
%   TransmitCorrelationMatrix  - Optional. Transmit spatial correlation
%                                matrix (or 3-D array) (default 1)
%   ReceiveCorrelationMatrix   - Optional. Receive spatial correlation
%                                matrix (or 3-D array) (default [1 0; 0 1])
%   TransmitPolarizationAngles - Optional. Transmit polarization slant
%                                angles in degrees (default [45 -45])
%   ReceivePolarizationAngles  - Optional. Receive polarization slant
%                                angles in degrees (default [90 0])
%   XPR                        - Optional. Cross polarization power ratio
%                                (dB) (default 10)
%   SpatialCorrelationMatrix   - Optional. Combined correlation matrix
%                                (or 3-D array) (default [1 0; 0 1])
%   NormalizePathGains         - Optional. Normalize path gains (logical)
%                                (default true)
%   InitialTime                - Optional. Start time of fading process (s)
%                                (default 0)
%   NormalizeChannelOutputs    - Optional. Normalize channel outputs by
%                                number of receive antennas (logical) (default true)
%   NumTimeSamples             - Optional. Number of time samples (default 30720)
%   NumSinusoids               - Optional. Number of sinusoids (default 48)
%   RandomStream               - Optional. Source of random number stream
%                                ("Global stream", "mt19937ar with seed" (default))
%   Seed                       - Optional. Initial seed of mt19937ar random
%                                number stream generator (default 73)
%
%   Note that the maximum Doppler shift caused by the movement of both
%   mobile and satellite, must be less than one-tenth of SampleRate.
%
%   % Example:
%   % Setup the NTN narrowband channel for an urban environment. Consider a
%   % LEO satellite operating in Ka-band at an altitude of 1500 km and
%   % speed of 7.1172 km/s. Assume a static UE and elevation angle of
%   % 50 degrees.
%
%   chanParams = struct;
%   chanParams.NTNChannelType = "Narrowband";
%   chanParams.Environment = "Urban";
%   chanParams.CarrierFrequency = 20e9;     % In Hz
%   chanParams.ElevationAngle = 50;         % In degrees
%   chanParams.SatelliteSpeed = 7117.2;     % In m/s
%   chanParams.SatelliteAltitude = 1500000; % In m
%   chanParams.MobileSpeed = 0;             % In m/s
%   chanParams.AzimuthOrientation = 0;      % In degrees
%
%   ntnNarrowbandChan = HelperSetupNTNChannel(chanParams);
%
%   See also HelperGenerateNTNChannel, nrTDLChannel, p681LMSChannel,
%   comm.ChannelFilter.

%   Copyright 2021-2023 The MathWorks, Inc.

    % Check the mandatory fields
    chanParams = validateMandatoryFields(chanParams);

    % Assign temporary variables for carrier frequency and maximum Doppler
    % shift due to mobile movement
    fc = double(chanParams.CarrierFrequency);
    maxDoppler = (double(chanParams.MobileSpeed)*fc)/physconst('LightSpeed');

    if strcmpi(chanParams.NTNChannelType,'TDL')
        % Construct nrTDLChannel
        baseChan = nrTDLChannel;
        baseChan = passign(chanParams,baseChan,{'DelaySpread'}, ...
            ~strcmp(chanParams.DelayProfile,'Custom'));
        baseChan.DelayProfile = 'Custom';
        baseChan = passign(chanParams,baseChan,{'FadingDistribution'}, ...
            strcmp(chanParams.DelayProfile,'Custom'));

        % Get the delay-profiles of the NTN channel and assign them to
        % nrTDLChannel properties
        hasLOSPath = (any(strcmpi(chanParams.DelayProfile,{'NTN-TDL-C','NTN-TDL-D'}))) ...
        || (strcmpi(baseChan.FadingDistribution,'Rician'));
        if strcmpi(chanParams.DelayProfile,'Custom')
            baseChan = passign(chanParams,baseChan,{'PathDelays'});
            baseChan = passign(chanParams,baseChan,{'AveragePathGains'});
            baseChan = passign(chanParams,baseChan,{'KFactorFirstTap'},hasLOSPath);
        else
            desiredKFactor = NaN;
            if hasLOSPath
                % Update to a LOS channel, to check the values of
                % KFactorScaling and KFactor properties
                baseChan.DelayProfile = 'TDL-D';
                baseChan = passign(chanParams,baseChan,{'KFactorScaling'});
                baseChan = passign(chanParams,baseChan,{'KFactor'}, ...
                    baseChan.KFactorScaling);
                if baseChan.KFactorScaling
                    desiredKFactor = double(baseChan.KFactor);
                else
                    desiredKFactor = NaN;
                end
                % Update to Custom profile
                baseChan.DelayProfile = 'Custom';
            end
            switch lower(chanParams.DelayProfile)
                case 'ntn-tdl-a'
                    pdp = [0 0; ...
                        1.0811 -4.675; ...
                        2.8416 -6.482];
                case 'ntn-tdl-b'
                    pdp = [0 0; ...
                        0.7249 -1.973; ...
                        0.7410 -4.332; ...
                        5.7392 -11.914];
                case 'ntn-tdl-c'
                    pdp = [0 -0.394; ...
                        0 -10.618; ...
                        14.8124 -23.373];
                otherwise
                    pdp = [0 -0.284; ...
                        0 -11.991; ...
                        0.5596 -9.887; ...
                        7.3340 -16.771];
            end
            % Perform delay and K-factor scaling
            pdp = double(wireless.internal.channelmodels.scaleDelaysAndKFactor( ...
                pdp,desiredKFactor,baseChan.DelaySpread));

            pathDelays = pdp(:,1).'; % 1st column is delay
            pathGains = pdp(:,2).';  % 2nd column is power

            if hasLOSPath
                % Remove 2nd entry of delay profile, corresponding to the
                % Rayleigh part of the Rician path. The Rician path is now
                % captured through first entry of path gains and the
                % K-factor of first tap
                baseChan.FadingDistribution = 'Rician';
                baseChan.KFactorFirstTap = pathGains(1) - pathGains(2);
                pathGainsOut = [10*log10(sum(10.^(pathGains(1:2)/10))) pathGains(3:end)];
                pathDelaysOut = pathDelays([1 3:end]);
            else
                baseChan.FadingDistribution = 'Rayleigh';
                pathGainsOut = pathGains;
                pathDelaysOut = pathDelays;
            end
            baseChan.PathDelays = pathDelaysOut;
            baseChan.AveragePathGains = pathGainsOut;
        end

        % Assign other properties of nrTDLChannel
        baseChan.ChannelFiltering = false; % Set ChannelFiltering to false
        baseChan = passign(chanParams,baseChan,{'OutputDataType'});
        baseChan.MaximumDopplerShift = maxDoppler;
        [baseChan,mimoFlag] = passign(chanParams,baseChan,{'MIMOCorrelation'});
        customMIMOProfile = mimoFlag && strcmp(chanParams.MIMOCorrelation,'Custom');
        [baseChan,polarizationFlag] = passign(chanParams,baseChan,{'Polarization'});
        customPolarization = polarizationFlag && strcmp(chanParams.Polarization,'Custom');
        customMIMOAndCrossPolar = ...
            customMIMOProfile && strcmp(baseChan.Polarization,'Cross-Polar');
        inFieldNames = {'TransmissionDirection', 'NumReceiveAntennas', ...
            'NumTransmitAntennas', 'TransmitCorrelationMatrix', 'ReceiveCorrelationMatrix', ...
            'TransmitPolarizationAngles', 'ReceivePolarizationAngles', 'XPR', ...
            'SpatialCorrelationMatrix', 'NormalizePathGains', 'NumSinusoids', ...
            'InitialTime', 'RandomStream', 'NormalizeChannelOutputs', ...
            'NumTimeSamples', 'SampleRate'};
        acCond = [repmat(~customMIMOProfile,1,2)  ~(customMIMOProfile && ~customPolarization) ...
            repmat(~(~customMIMOProfile || customPolarization),1,2) ...
            repmat(customMIMOAndCrossPolar,1,3) (customMIMOProfile && customPolarization) ...
            ones(1,7)];
        baseChan = passign(chanParams,baseChan,inFieldNames,acCond);
        baseChan = passign(chanParams,baseChan,{'Seed'}, ...
            strcmp(baseChan.RandomStream,'mt19937ar with seed'));
        pathDelays = baseChan.PathDelays;
        chanName = "NTN TDL with " + chanParams.DelayProfile + " delay profile";
    else
        % Construct p681LMSChannel
        baseChan = p681LMSChannel;
        baseChan.CarrierFrequency = chanParams.CarrierFrequency;
        baseChan.MobileSpeed = chanParams.MobileSpeed;
        baseChan.ElevationAngle = chanParams.ElevationAngle;
        baseChan.AzimuthOrientation = chanParams.AzimuthOrientation;
        baseChan.Environment = chanParams.Environment;
        % Parse the optional fields
        baseChan.ChannelFiltering = false; % Set ChannelFiltering to false
        baseChan = passign(chanParams,baseChan,{'OutputDataType'});
        inFieldNames = {'StateDistribution', 'MinStateDuration', ...
            'DirectPathDistribution', 'MultipathPowerCoefficients', ...
            'StandardDeviationCoefficients', 'DirectPathCorrelationDistance', ...
            'TransitionLengthCoefficients', 'StateProbabilityRange', 'NumSamples', ...
            'FadingTechnique', 'RandomStream', 'SampleRate', 'InitialState'};
        customEnvi = strcmpi(baseChan.Environment,'Custom');
        acCheck = [repmat(customEnvi,1,8) ones(1,5)];
        baseChan = passign(chanParams,baseChan,inFieldNames,acCheck);
        baseChan = passign(chanParams,baseChan,{'NumSinusoids'}, ...
            strcmpi(baseChan.FadingTechnique,'sum of sinusoids'));
        baseChan = passign(chanParams,baseChan,{'Seed'}, ...
            strcmpi(baseChan.RandomStream,'mt19937ar with seed'));
        pathDelays = 0;
        chanName = "NTN narrowband with " + baseChan.Environment + " environment";
    end

    % Calculate the Doppler shift due to satellite movement
    fdSat = satcom.internal.dopplerShift(fc,double(chanParams.SatelliteSpeed), ...
        double(chanParams.ElevationAngle),double(chanParams.SatelliteAltitude));
    % Check the maximum Doppler shift and sample rate
    if ((maxDoppler+abs(fdSat)) >= (baseChan.SampleRate/10))
        error("satcom:HelperSetupNTNChannel:MaxDoppler",...
            "The maximum Doppler shift (%d) due to mobile and satellite "+ ...
            "movement, must be less than one-tenth of SampleRate.", ...
            (maxDoppler+abs(fdSat)));
    end

    % Set the channel filter
    chanFilt = comm.ChannelFilter(...
                'SampleRate',baseChan.SampleRate,'PathDelays',pathDelays, ...
                'NormalizeChannelOutputs',false);

    % Set the output structure
    ntnChan = struct;
    ntnChan.ChannelName = chanName;
    ntnChan.BaseChannel = baseChan;
    ntnChan.ChannelFilter = chanFilt;
    ntnChan.SatelliteDopplerShift = fdSat;
    ntnChan.MobileDopplerSpread = maxDoppler;

end

function chanParams = validateMandatoryFields(chanParams)
% Validate the mandatory fields

    if isfield(chanParams,'NTNChannelType') && ...
            isfield(chanParams,'CarrierFrequency') && ...
            isfield(chanParams,'SatelliteSpeed') && ...
            isfield(chanParams,'SatelliteAltitude') && ...
            isfield(chanParams,'MobileSpeed') && ...
            isfield(chanParams,'ElevationAngle')
        chanParams.NTNChannelType = validatestring(chanParams.NTNChannelType, ...
            {'Narrowband','TDL'},'','CHANPARAMS.NTNChannelType');
        validateattributes(chanParams.CarrierFrequency, {'double'}, ...
                {'scalar','real','finite','nonnegative'},'','CHANPARAMS.CarrierFrequency');
        validateattributes(chanParams.SatelliteSpeed, {'double'}, ...
                {'scalar','real','finite'},'','CHANPARAMS.SatelliteSpeed');
        validateattributes(chanParams.SatelliteAltitude, {'double'}, ...
                {'scalar','real','finite','nonnegative'},'','CHANPARAMS.SatelliteAltitude');
        validateattributes(chanParams.MobileSpeed, {'double'}, ...
                {'scalar','real','finite','nonnegative'},'','CHANPARAMS.MobileSpeed');
        validateattributes(chanParams.ElevationAngle, {'double'}, ...
                {'scalar','real','finite'},'','CHANPARAMS.ElevationAngle');
        if strcmpi(chanParams.NTNChannelType,'TDL')
            if isfield(chanParams,'DelayProfile')
                chanParams.DelayProfile = validatestring(chanParams.DelayProfile, ...
                    {'NTN-TDL-A','NTN-TDL-B','NTN-TDL-C','NTN-TDL-D','Custom'}, ...
                    '','CHANPARAMS.DelayProfile');
            else
                error("satcom:HelperSetupNTNChannel:TDLRequiredFields", ...
                    "CHANPARAMS must have DelayProfile field when " + ...
                    "NTNChannelType field is set to TDL.");
            end
        else % Narrowband
            if isfield(chanParams,'AzimuthOrientation') && ...
                    isfield(chanParams,'Environment')
                validateattributes(chanParams.AzimuthOrientation, {'double'}, ...
                {'scalar','real','finite'},'','CHANPARAMS.AzimuthOrientation');
                chanParams.Environment = validatestring(chanParams.Environment, ...
                    {'Urban', 'Suburban', 'RuralWooded', 'Village', 'Train', ...
                    'Residential', 'Highway', 'Rural', 'Custom'},'', ...
                    'CHANPARAMS.Environment');
            else
                error("satcom:HelperSetupNTNChannel:NarrowbandRequiredFields", ...
                    "CHANPARAMS must have AzimuthOrientation and " + ...
                    "Environment fields, when NTNChannelType field is set to Narrowband.");
            end
        end
    else
        error("satcom:HelperSetupNTNChannel:CommonRequiredFields", ...
            "CHANPARAMS must have the mandatory fields: NTNChannelType, " ...
            + "CarrierFrequency, SatelliteSpeed, SatelliteAltitude, " ...
            + "ElevationAngle, and MobileSpeed.");
    end

end

function [o,cond] = passign(s,o,f,ac)
% Assign the fields of structure to the appropriate object properties

    cond = isfield(s,f);
    if nargin == 4
        cond = cond & ac;
    end

    for i = 1:length(cond)
        if cond(i)
            o.(f{i}) = s.(f{i});
        end
    end

end
