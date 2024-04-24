%% functions for transform antenna power from (dBm) to (W)

% W => dBm

function [UE_antenna_dBm] = mW_to_dBm(UE_antenna_power)

    % Transmit power : 200mW (paper parameters setup)
    % 200mW = 0.2W
    % dBm = 10 * log(P/1mW)
    % 23 dBm = 10 * log10(200)

    UE_antenna_dBm = 10*log(UE_antenna_power);

end