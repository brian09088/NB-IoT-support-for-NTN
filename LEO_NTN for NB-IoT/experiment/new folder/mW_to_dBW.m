%% functions for transform antenna power from (mW) to (dBW)

% mW => dBW

function [dBW] = mW_to_dBW(txpower)

    % Transmit power : 200mW (paper parameters setup)
    % 200mW = 0.2W
    % dBm = 10 * log(P/1mW)
    % 23 dBm = 10 * log10(200)
    % dBm = dBW + 30
    % dBW = -7

    dBW = 10*log(txpower) - 30;

end