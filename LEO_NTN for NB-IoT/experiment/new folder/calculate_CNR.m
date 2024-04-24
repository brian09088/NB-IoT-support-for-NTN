% function for calculate CNR

function [cn] = calculate_CNR(txpower_dBW ,fc_GHz , Transmit_antenna_gain, d, G_T, data_rate)

    cfg = satelliteCNRConfig;
    cfg.TransmitterPower = txpower_dBW;                         % in dBW ()
    cfg.TransmitterSystemLoss = 0;                              % in dB (FSPL)
    cfg.TransmitterAntennaGain = Transmit_antenna_gain;         % in dBi
    cfg.Distance = d/1000;                                      % in km (LEO sat altitude)
    cfg.Frequency = fc_GHz;                                     % in GHz
    
    % Here, miscellaneous losses include polarization loss, interference 
    % loss, and antenna mispointing loss, respectively.
    polLoss = 0;
    intLoss = 0;
    antLoss = 0;
    cfg.MiscellaneousLoss =  polLoss + intLoss + antLoss; % in dB
    cfg.GainToNoiseTemperatureRatio = G_T;                % in dB/K
    cfg.ReceiverSystemLoss = 0;                           % in dB
    % Among current NTN target performance scenarios
    % IoT connectivity requires a data rate of 10 kbps or more
    cfg.BitRate = data_rate;                                   % in Mbps
    
    disp(cfg)
    
    [cn,info] = satelliteCNR(cfg)

end


    

