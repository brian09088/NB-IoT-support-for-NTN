% function for calculate doppler compensation

function [dopp_range] = dopp_compensation(FFT_size, theta, L)

    dopp_range = (FFT_size * theta) / (2 * pi * L);

end