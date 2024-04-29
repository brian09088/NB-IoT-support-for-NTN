function chanOut = setupNTNChannel(channel,sampleRate)
% Setup NTN channel

    % Assign temporary variables for carrier frequency and maximum Doppler
    % shift due to mobile movement
    fc = double(channel.CarrierFrequency);
    c = physconst("LightSpeed");
    maxDoppler = (double(channel.MobileSpeed)*fc)/c;
    elevAngle = double(channel.ElevationAngle);
    h = double(channel.SatelliteAltitude);
    v = double(channel.SatelliteSpeed);
    % Calculate the Doppler shift due to satellite movement
    maxDopplerSat = satcom.internal.dopplerShift(fc,v,elevAngle,h);
    % Check the maximum Doppler shift and sample rate
    if ((maxDoppler+maxDopplerSat) >= (sampleRate/10))
        error("satcom:setupNTNChannel:MaxDoppler", ...
            "The maximum Doppler shift (%d Hz) due to mobile and satellite " + ...
            "movement, must be less than %d Hz which is one-tenth of SampleRate.", ...
            (maxDoppler + maxDopplerSat),sampleRate/10)
    end

    chanOut = struct;
    chanTypeLower = lower(channel.NTNChannelType);
    if chanTypeLower == "etsi rician"
        channelName = "ETSI Rician";
        baseChannel = etsiRicianChannel;
        baseChannel.SampleRate = sampleRate;
        baseChannel.KFactor = channel.KFactor;
        baseChannel.MaximumDopplerShift = maxDoppler;
    elseif chanTypeLower == "itu-r p.681"
        channelName = "ITU-R P.681";
        baseChannel = p681LMSChannel;
        baseChannel.SampleRate = sampleRate;
        baseChannel.Environment = channel.Environment;
        baseChannel.CarrierFrequency = channel.CarrierFrequency;
        baseChannel.MobileSpeed = channel.MobileSpeed;
        baseChannel.ElevationAngle = channel.ElevationAngle;
        baseChannel.AzimuthOrientation = channel.AzimuthOrientation;
        baseChannel.FadingTechnique = "Sum of sinusoids";
    end
    baseChannel.RandomStream = "mt19937ar with seed";
    baseChannel.Seed = channel.Seed;

    % Set the channel filter
    chanFilt = comm.ChannelFilter( ...
                SampleRate=sampleRate,PathDelays=0, ...
                NormalizeChannelOutputs=false);

    % Set the output structure
    chanOut.ChannelName = channelName;
    chanOut.CarrierFrequency = fc;
    chanOut.SatelliteSpeed = v;
    chanOut.SatelliteAltitude = h;
    chanOut.ElevationAngle = elevAngle;
    chanOut.BaseChannel = baseChannel;
    chanOut.SatelliteDopplerShift = maxDopplerSat;
    chanOut.ChannelFilter = chanFilt;

end

function [out,sampleTimes] = generateNTNChannel(channel,in)
% Generate NTN channel

    % Get the channel information before channel processing
    prevInfo = info(channel.BaseChannel);
    numSamplesStart = prevInfo.NumSamplesProcessed;

    % Get the path gains of base channel
    [~,pathGainsBase] = channel.BaseChannel(in);

    % Get the channel information after channel processing
    postInfo = info(channel.BaseChannel);
    numSamplesEnd = postInfo.NumSamplesProcessed;

    % Get the channel sample times
    sampleTimes = (numSamplesStart:(numSamplesEnd-1)).'/channel.BaseChannel.SampleRate;

    % Apply satellite Doppler shift to the base channel path gains
    pathGains = pathGainsBase.*exp(1i*2*pi*channel.SatelliteDopplerShift*sampleTimes);

    % Perform channel filtering
    out = channel.ChannelFilter(in,pathGains);

end

function out = compensateDopplerShift(enb,inWave,foffset,flag)
% Perform Doppler shift correction

    if flag
        % Correct frequency offset
        out = lteFrequencyCorrect(enb,inWave,foffset);
    else
        out = inWave;
    end

end

function out = estimateDopplerShiftUsingRS(enb,rxWave,refInd, ...
    refSym,flag)
% Estimate the Doppler shift using NPSS

    if flag
        % Set the Windowing field to 0, as this information is not known at
        % the receiver
        enb.Windowing = 0;
        ofdmInfo = lteSCFDMAInfo(enb);
        K = 12;                             % Number of subcarriers     
        L = 14;                             % Number of OFDM symbols in slot

        % Initialize temporary variables
        rxWave1 = [rxWave; zeros((mod(size(rxWave,1),2)),1)]; % Append zero, if required
        rxLen = size(rxWave1,1);

        % Generate reference waveform
        refGrid = complex(zeros([K L]));
        refGrid(refInd) = refSym;
        refWave = lteSCFDMAModulate(enb,refGrid);
        refWave = [refWave; zeros((rxLen-size(refWave,1)),1)];

        % Compute the correlation of received waveform with reference
        % waveform
        x_wave = rxWave1.*conj(refWave);

        % Compute FFT of the resultant waveform
        x_fft = fftshift(fft(x_wave));

        % FFT bin values
        fftBinValues = (-rxLen/2:(rxLen/2-1))*(ofdmInfo.SamplingRate/rxLen);

        % Use the FFT bin index corresponding to the maximum FFT value.
        % The FFT bin value corresponding to this bin index is the integer
        % frequency offset.
        [~,binIndex] = max(x_fft);
        out = fftBinValues(binIndex);
    else
        out = 0;
    end

end

function varargout = paMemorylessNonlinearity(paModel,varargin)
% Apply power amplifier nonlinearity (TR 38.803)
% out = paMemorylessNonlinearity(paModel,in,enable) returns the
% impaired output.
% paMemorylessNonlinearity(paModel) returns the plot with the gain and
% phase characteristics of the power amplifier

    if nargin == 1
        in_NoScale = randn(1e6,1)+1j*randn(1e6,1);
        scaleFactor = 1/sqrt(2);
        enable = 1;
    else
        in_NoScale = varargin{1};
        scaleFactor = 1;
        enable = varargin{2};
    end

    if enable
        in = scaleFactor*in_NoScale;
        if isa(paModel,"comm.MemorylessNonlinearity")
            % paModel is a comm.MemorylessNonlinearity System object
            out = paModel(in);
            paModelName = "";
        else
            absIn = abs(in);
            paModelName = paModel;
            switch lower(paModel)
                case "2.1ghz gaas"
                    % 2.1GHz GaAs
                    out = (-0.618347-0.785905i) * in + (2.0831-1.69506i) * in .* absIn.^(2) + ...
                        (-14.7229+16.8335i) * in .* absIn.^(2*2) + (61.6423-76.9171i) * in .* absIn.^(2*3) + ...
                        (-145.139+184.765i) * in .* absIn.^(2*4) + (190.61-239.371i)* in .* absIn.^(2*5) + ...
                        (-130.184+158.957i) * in .* absIn.^(2*6) + (36.0047-42.5192i) * in .* absIn.^(2*7);
                otherwise
                    % 2.1GHz GaN
                    out = (0.999952-0.00981788i) * in + (-0.0618171+0.118845i) * in .* absIn.^(2) + ...
                        (-1.69917-0.464933i) * in .* absIn.^(2*2) + (3.27962+0.829737i) * in .* absIn.^(2*3) + ...
                        (-1.80821-0.454331i) * in .* absIn.^(2*4);
            end
        end
    else
        out = in_NoScale;
    end

    if nargout > 0
        varargout{1} = out;
    end

    if nargin == 1 || (nargout == 0)
        % Gain Plot
        inpPower = 20*log10(absIn);
        gain = 20*log10(abs(out))-inpPower;
        figure
        subplot(211)
        plot(inpPower,gain,".")
        grid on
        ylim([-Inf 1])
        xlim([-30 0])
        xlabel("Normalized input power (dB)")
        ylabel("Gain (dB)")
        title("Gain Characteristics of PA Model " + paModelName)

        % Phase Plot
        phase = angle(out.*conj(in))*180/pi;
        subplot(212)
        plot(inpPower,phase,".")
        grid on
        xlim([-30 0])
        xlabel("Normalized input power (dB)")
        ylabel("Phase (deg)")
        title("Phase Characteristics of PA Model " + paModelName)
    end

end

function paChar = getDefaultCustomPA()
% The operating specifications for the LDMOS-based Doherty amplifier are:
% * A frequency of 2110 MHz
% * A peak power of 300 W
% * A small signal gain of 61 dB
% Each row in HAV08_Table specifies Pin (dBm), gain (dB), and phase shift
% (degrees) as derived from figure 4 of Hammi, Oualid, et al. "Power
% amplifiers' model assessment and memory effects intensity quantification
% using memoryless post-compensation technique." IEEE Transactions on
% Microwave Theory and Techniques 56.12 (2008): 3170-3179.

    HAV08_Table = ...
        [-35,60.53,0.01;
        -34,60.53,0.01;
        -33,60.53,0.08;
        -32,60.54,0.08;
        -31,60.55,0.1;
        -30,60.56,0.08;
        -29,60.57,0.14;
        -28,60.59,0.19;
        -27,60.6,0.23;
        -26,60.64,0.21;
        -25,60.69,0.28;
        -24,60.76,0.21;
        -23,60.85,0.12;
        -22,60.97,0.08;
        -21,61.12,-0.13;
        -20,61.31,-0.44;
        -19,61.52,-0.94;
        -18,61.76,-1.59;
        -17,62.01,-2.73;
        -16,62.25,-4.31;
        -15,62.47,-6.85;
        -14,62.56,-9.82;
        -13,62.47,-12.29;
        -12,62.31,-13.82;
        -11,62.2,-15.03;
        -10,62.15,-16.27;
        -9,62,-18.05;
        -8,61.53,-20.21;
        -7,60.93,-23.38;
        -6,60.2,-26.64;
        -5,59.38,-28.75];
    % Convert the second column of the HAV08_Table from gain to Pout for
    % use by the memoryless nonlinearity System object.
    paChar = HAV08_Table;
    paChar(:,2) = paChar(:,1) + paChar(:,2);

end