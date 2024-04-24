%% function for calculate contact angle

function [phi] = calculate_phi(gs_lat, gs_lon)
    
    % (1) get satellite location from tle file (starlink 1007~1009),...
    % calculate with timestamp: starlink-1007 launch date(2019/11/11)
    % (2) then ccalculate phi from lat lon
    % (3) return phi for following calculate_d
    
    % defined satellite constellation
    % space-X starlink satellite dconstellation
    % ex: starlink_tle.csv, starlink_1007.csv
    sat_csv = uigetfile('*.csv','Choose a .csv file');
    
    % 將TLE透過satellite toolbox導入
    location = satelliteTLE(sat_csv);
    
    sat_lat = location(:, 1);
    sat_lon = location(:, 2);
    sat_alt = location(:, 3);

    phi = acosd(cos(sat_lat)*cos(gs_lat)*cos(sat_lon-gs_lon) + sin(sat_lat)*sin(sat_lon)); 
    % (degree)

end