%% function for power transfer 

% from dBm to mW

function [mW] = dBW_to_mW(dBW)

    dBm = dBW - 30;
    mW = 10 ^ (dBm/10);

end
