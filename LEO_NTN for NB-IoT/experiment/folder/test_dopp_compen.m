FFT_size = 128;    % Fast Fourier Transform size (n=128)
n = 10;            % number of symbols between the DMRS symbols
Lcp = num2cell([10, repmat(9, 1, n-1)]);  % Assign 10 to the first symbol and 9 to the rest

% Number of time samples between channel estimation symbols
% L = calculate_L(Lcp, n, FFT_size);   
L = 960;

% Phase difference between channel estimation symbols, theta = 0
theta = calculate_theta();  

% range of the Doppler shift compensation
dopp_range = dopp_compensation(FFT_size, theta, L);

% maximum phase difference between DMRS symbols = pi, dopp_range = 1/15
% SCS = 15 kHz, fd_max = 1kHz
SCS = 15;
fd_max = SCS * dopp_range;

fprintf("fd_max = %d kHz \n", fd_max);

function theta = calculate_theta()
    % phase difference between channel estimation symbols
    theta = pi;
end

function [L] = calculate_L(Lcp, n, FFT_size)
    % Initialize variable to store the total number of time samples
    total_Lcp = 0;
    
    % Loop through each symbol
    for i = 1:n
        % Accumulate the cyclic prefix time samples
        total_Lcp = total_Lcp + Lcp{i};
    end
    
    % Calculate the total number of time samples between symbols
    L = total_Lcp + n * FFT_size;
end

function [dopp_range] = dopp_compensation(FFT_size, theta, L)
    dopp_range = (FFT_size * theta) / (2 * pi * L);
end
