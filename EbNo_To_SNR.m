% This function converts the Eb/No to SNR in dB

function [SNR] = EbNo_To_SNR(EbNo, BitRate , Bandwidth)

% EbNo is in dB
% BitRate and Bandwidth must have the same units

SNR = EbNo + 10*log10(BitRate/Bandwidth); % SNR in dB

end
