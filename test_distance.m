%% for function test

gs_name = 'south-korea';
gs_lat = 36.35080;
gs_lon = 127.30122;
gs_location = [gs_lat, gs_lon];

Re = 6371;          % Radius of Earth(m)
h0 = 650;           % altitude of satellite(m)
alpha = 53.1;            % 發射傾角 degrees (>10 degree)

phi = calculate_phi(gs_lat, gs_lon);
% distance between sat & UE (meters)
d = cal_distance(Re, phi, h0);
fprintf("d = %.3d\n",d)

% distance between sat & UE (meters)
d1 = calculate_d(Re, alpha, h0);  
d2 = calculate_d(Re, alpha, h0);  

% fprintf("use paper1 formula : d = %.3d\n",d1);
% fprintf("use paper2 formula : d = %.3d\n",d1);
fprintf("use paper 1 formula : d = %.3d\n",d1)
fprintf("use paper 2 formula : d = %.3d\n",d2)


function [lat,lon,alt] = satelliteTLE(sat_csv)
    
    startTime = datetime(2021,11,11,8,30,0);
    sampletime = 5736;              % seconds
    stopTime = startTime + sampletime;      
    sc = satelliteScenario(startTime,stopTime,sampletime);

    sat = satellite(sc,sat_csv);
    time = datetime(2021,11,11,9,30,0);
    [lat,lon,alt] = states(sat,time,"CoordinateFrame","geographic");

end

function [d1,d2] = calculate_d(Re, alpha, h0)

    d1 = sqrt(Re^2 * (sin(alpha)^2) + h0^2 + 2*h0*Re) - Re * sin(alpha);
    d2 = sqrt(Re^2 + (Re+h0)^2 - 2*Re*(Re+h0)*cos(alpha-90));

end

% function for calculate contact angle

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
    
    sat_lat = location(1);
    sat_lon = location(2);

    phi = acosd(cos(sat_lat)*cos(gs_lat)*cos(sat_lon-gs_lon) + sin(sat_lat)*sin(sat_lon)); 
    % (degree)

end