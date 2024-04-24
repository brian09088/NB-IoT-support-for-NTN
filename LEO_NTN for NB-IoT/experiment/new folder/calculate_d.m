%% function use for calculate distance between satellite and UE

% use inclination(elevation angle) to calculate

function [d] = calculate_d(Re, alpha, h0)

    %% according paper1(Performance analysis of NB-IoT...) formula
    d = sqrt(Re^2 * (sin(alpha)^2) + (h0)^2 + 2*h0*Re) - Re*sin(alpha);

    %% according paper2(Link budget analysis for satellite_based NB-IoT...) formula
    % ✅ test output same as paper1 
    % d = -Re * sin(alpha) + sqrt(Re^2*sin(alpha)^2 + h0 + 2*Re*h0);
    
    %% use law of cosine to calculate
    % d = sqrt(Re^2 + (Re+h0)^2 - 2*Re*(Re+h0)*cos(alpha));
    % d = sqrt(Re^2 + (Re+h0)^2 - 2*Re*(Re+h0)*cos(alpha-90));
    % d = sqrt(Re^2 + (Re+h0)^2 - 2*Re*(Re+h0)*cos(90-alpha));

end