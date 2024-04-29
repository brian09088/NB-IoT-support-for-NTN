%% test adding quoted file path

% sample output
% >> test
% E:\MATLAB\norad.tle
% 'E:\MATLAB\norad.tle'

% Your TLE file path
tleFilePath = 'E:\MATLAB\norad.tle';

% Enclose the path in single quotes
quotedTleFilePath = ['''', tleFilePath, ''''];

% Display the quoted path to confirm
disp(tleFilePath);
disp(quotedTleFilePath);
