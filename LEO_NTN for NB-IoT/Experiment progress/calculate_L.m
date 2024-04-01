% function for calculate L

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