cfg = satelliteCNRConfig;
cfg.TransmitterPower = -7;                         % in dBW ()
cfg.TransmitterSystemLoss = 0;                              % in dB (FSPL)
cfg.TransmitterAntennaGain = 0;         % in dBi
cfg.Distance = 1397;                                         % in km (LEO sat altitude)
cfg.Frequency = 2;                                     % in GHz

% Here, miscellaneous losses include polarization loss, interference 
% loss, and antenna mispointing loss, respectively.
polLoss = 0;
intLoss = 0;
antLoss = 0;
cfg.MiscellaneousLoss =  polLoss + intLoss + antLoss; % in dB
cfg.GainToNoiseTemperatureRatio = 1.1;                % in dB/K
cfg.ReceiverSystemLoss = 0;                           % in dB
% Among current NTN target performance scenarios
% IoT connectivity requires a data rate of 10 kbps or more
cfg.BitRate = 0.01;                                   % in Mbps

disp(cfg)

[cn,info] = satelliteCNR(cfg)