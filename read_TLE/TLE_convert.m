clc,clear,close all

%% convert TLE file to orbit parameters

% get constellation TLE file from csv/txt... download from websites
% sat_tle = uigetfile('*.csv','Choose a .csv file');
% sat_tle = tleread('***.tle');

fullPath = mfilename('fullpath');   % % Get full path of current function
[currentFolder, ~, ~] = fileparts(fullPath);   % Get the folder containing function
tleFileName = 'norad.tle';  % Define name of TLE file
tleFilePath = fullfile(currentFolder, tleFileName); % Construct full path

% convert from TLE files, get 6 orbit parameter elements
convert_TLE(tleFilePath)

% Create satellite in MATLAB
sat = satellite(sc, semiMajorAxis, eccentricity, inclination, ...
    rightAscensionOfAscendingNode, argumentOfPeriapsis, ...
    trueAnomaly, Name="NORAD_1 ~ NORAD_66");

show(sat)

function [mu, semi_major_axis, eccentricity, ...
    inclination, RA_of_asc_node, Arg_of_perigee, ...
    Mean_anomaly, Mean_motion] = convert_TLE(tleFilePath)

    % Earth's gravitational parameter
    % mu = 3.986e14;      % m^3/s^2
    
    % semi_major_axis;    % meters
    % eccentricity;
    % inclination;        % deg
    % RA_of_asc_node;     % deg
    % Arg_of_perigee;     % deg
    % Mean_anomaly;       % deg
    % Mean_motion;        % rev/day

    % wait me keep doing this part
    % convert formula...
    
    % Your TLE file path
    tleFilePath = 'E:\MATLAB\碩士論文\Brian_Su\read_TLE\norad.tle';

    % Enclose the path in single quotes
    quotedTleFilePath = ['''', tleFilePath, ''''];
    
    readtle(quotedTleFilePath)

end