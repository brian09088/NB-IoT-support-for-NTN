%% function use for calculate distance between satellite and UE

% use central angle : phi

function [d] = cal_distance(Re, phi, h0)

    d = h0*(1+(Re/h0)^2 - 2*(Re/h0)*cos(phi))^(1/2);

end