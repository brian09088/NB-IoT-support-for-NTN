%% functions for transform antenna power from (dBm) to (W)

% W => dBm

function [UE_antenna_dBm] = mW_to_dBm(UE_antenna_power)

    % Transmit power : 200mW (paper parameters setup)
    % 200mW = 0.2W
    % dBm = 10 * log(W)
    % -6.99 dBm = 10 * log(0.2)

    UE_antenna_dBm = 10*log(UE_antenna_power);

end