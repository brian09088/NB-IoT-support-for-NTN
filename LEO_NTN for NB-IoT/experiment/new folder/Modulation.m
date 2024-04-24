% This function computes the BER for different types of mod_orders. 
% For a given mod_order, mod_order order, and BER specification, it gives the required Eb/NO.

function [EbNo_spec] = Modulation(mod_order , ber_spec)

EbNo = 0:0.01:25;

switch mod_order
                    
    case {"OQPSK" ,"oqpsk"}
        dataenc = input('Differential encode (diff) or non differential (nondiff): ','s');
        BER = berawgn(EbNo, 'oqpsk' ,dataenc) ;

    case {"QPSK" ,"qpsk"}
    %   dataenc = input('Differential encode (diff) or non differential (nondiff): ','s');
    %   BER = berawgn(EbNo, 'qpsk' ,dataenc) ;

end

EbNo_spec = interp1(BER,EbNo, ber_spec) ;
semilogy (EbNo, BER)
hold on
semilogy (EbNo_spec,ber_spec, 'bp')
legend('BER', 'BER specification')
title ('BER vs Eb/No')
xlabel('Eb/No')
ylabel('BER')

end
