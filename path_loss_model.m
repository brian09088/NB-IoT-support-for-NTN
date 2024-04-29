%% function for calculate path loss

% PL = path loss(dB)
% FSPL = free space path loss (dB)
% APL = atmospheric loss 0.07 (dB)
% SMPL = shadowing margin loss 3.00 (dB)
% SLPL = scintillation loss 2.20 (dB)
% Polarization loss 0.00 dB
% ADPL = additional loss 0.00 (dB)
% PL = FSPL + APL + SMPL + SLPL + ADPL

function [PL] = path_loss_model(FSPL)

    APL = 0.07;
    SMPL = 3.00;
    SLPL = 2.20;
    ADPL = 0.00;
    PL = FSPL + APL + SMPL + SLPL + ADPL;

end