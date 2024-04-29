function varargout = HelperGenerateNTNChannel(ntnChan,in)
%HelperGenerateNTNChannel Generate NTN channel
%
%   Note: This is a helper and its API and/or functionality may change
%   in subsequent releases.
%
%   [PATHGAINS,SAMPLETIMES] = HelperGenerateNTNChannel(NTNCHAN) generates
%   the path gains, PATHGAINS, and sample times, SAMPLETIMES, of New Radio
%   (NR) non-terrestrial network (NTN) channel provided the input structure
%   NTNCHAN. NTNCHAN is a structure with fields:
%   BaseChannel           - Base channel to generate the path gains of
%                           channel with ChannelFiltering set to 0. Base
%                           channel is either p681LMSChannel or
%                           nrTDLChannel System object. In case of
%                           p681LMSChannel, the number of channel
%                           samples is defined by NumSamples property. In
%                           case of nrTDLChannel, the number of samples is
%                           defined by NumTimeSamples property. The NTN
%                           channel is configured to NTN narrowband
%                           channel, when base channel is p681LMSChannel.
%                           The NTN channel is configured to NTN TDL
%                           channel, when base channel is nrTDLChannel.
%   SatelliteDopplerShift - A scalar value indicating the doppler shift due
%                           to satellite movement.
%
%   [OUT,PATHGAINS,SAMPLETIMES] = HelperGenerateNTNChannel(NTNCHAN,IN)
%   filters the input signal, IN, through a NTN narrowband or TDL fading
%   channel and returns the result in OUT, given the input NTN channel
%   structure NTNCHAN. The number of channel samples is defined by the
%   first column of input signal X. NTNCHAN is a structure with fields:
%   BaseChannel           - Base channel to generate the path gains of
%                           channel with ChannelFiltering set to 0. Base
%                           channel is either p681LMSChannel or
%                           nrTDLChannel System object.
%   ChannelFilter         - Filter signal using multipath gains at
%                           specified path delays using this
%                           comm.ChannelFilter System object.
%   SatelliteDopplerShift - A scalar value indicating the doppler shift due
%                           to satellite movement.
%   For NTN TDL channel, the number of columns in input X must be same as
%   the number of transmit antennas that is configured in the BaseChannel
%   field of NTNCHAN. For NTN narrowband channel, number of columns in
%   input X must be 1.
%
%   % Example:
%   % Generate the NTN TDL channel for NTN-TDL-B delay profile. Consider a
%   % LEO satellite operating in Ka-band at an altitude of 1500 km and
%   % speed of 7.1172 km/s. Assume a static UE and elevation angle of 50
%   % degrees.
%
%   chanParams = struct;
%   chanParams.NTNChannelType = "TDL";
%   chanParams.DelayProfile = "NTN-TDL-B";
%   chanParams.CarrierFrequency = 20e9;     % In Hz
%   chanParams.ElevationAngle = 50;         % In degrees
%   chanParams.SatelliteSpeed = 7117.2;     % In m/s
%   chanParams.SatelliteAltitude = 1500000; % In m
%   chanParams.MobileSpeed = 0;             % In m/s
%
%   % Setup the NTN TDL channel
%   ntnTDLChan = HelperSetupNTNChannel(chanParams);
%
%   % Generate the path gains and sample times of NTN TDL channel
%   [pathGains,sampleTimes] = HelperGenerateNTNChannel(ntnTDLChan);
%
%   See also HelperSetupNTNChannel, nrTDLChannel, p681LMSChannel,
%   comm.ChannelFilter.

%   Copyright 2021-2023 The MathWorks, Inc.

    narginchk(1,2);
    if nargin > 1
        narginchk(2,2);
    end
    nargoutchk(0,nargin+1);

    % Get the path gains from base channel input
    if nargin == 2
        % Overwrite with the first dimension of input signal
        if isprop(ntnChan.BaseChannel,'NumTimeSamples')
            ntnChan.BaseChannel.NumTimeSamples = size(in,1);
        else
            ntnChan.BaseChannel.NumSamples = size(in,1);
        end
    end
    [pathGains,sampleTimes] = ntnChan.BaseChannel();

    % Apply doppler shift due to satellite motion for the path gains
    % generated from base channel
    pathGainsOut = ...
        pathGains.*exp(1i*2*pi*ntnChan.SatelliteDopplerShift(1)*sampleTimes);

    % Apply channel filtering
    if nargin == 2
        validateattributes(in,{'double','single'}, ...
            {'size',[nan size(pathGainsOut,3)]},'','IN');
        out = ntnChan.ChannelFilter(in,pathGainsOut);
        varargout = {out, pathGainsOut, sampleTimes};
    else
        varargout = {pathGainsOut, sampleTimes};
    end

end