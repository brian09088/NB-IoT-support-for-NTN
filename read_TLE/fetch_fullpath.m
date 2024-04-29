%% test for fetch TLE file under same folder of currently execution file

%% sample output
% >> test
% E:\MATLAB\碩士論文\Brian_Su\read_TLE\norad.tle

% Get the full path of the current script or function
fullPath = mfilename('fullpath');

% Get the folder containing the script/function
[currentFolder, ~, ~] = fileparts(fullPath);

% Define the name of the TLE file
tleFileName = 'norad.tle';

% Construct the full path to the TLE file
tleFilePath = fullfile(currentFolder, tleFileName);