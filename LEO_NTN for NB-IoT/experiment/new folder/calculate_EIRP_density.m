%% function for calculate EIRP density

function [EIRP_density] = calculate_EIRP_density(EIRP, bandwidth_kHz)

    % EIRP_density (dBW/MHz) = EIRP(dBW) - 10*log10(Bandwidth in MHz)
    bandwidth_MHz = bandwidth_kHz / 1e3;
    EIRP_density = EIRP - 10*log10(bandwidth_MHz);

end
