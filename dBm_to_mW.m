%% function for power transfer 

% from dBm to mW

function [mW] = dBm_to_mW(dBm)

    mW = 10 ^ (dBm/10);

end
