function [ x_upsample ] = bits2symbols( packets, fsfd, m )
%   Transfer the data bits to symbols and upsample them to get the
%   transmitted signal
% Constellation or bit to symbol mapping
s = [exp(0) exp(1j*pi/4) exp(1j*3*pi/4) exp(1j*pi/2) exp(1j*7*pi/4) exp(1j*3*pi/2) exp(1j*pi) exp(1j*5*pi/4)]; % Constellation 1 - QPSK/4-QAM
                                                        % s = exp(1i*((0:3)*pi/2 + pi/4)); % Constellation 1 - same constellation generated as PSK
% scatterplot(s); grid on;                            % Constellation visualization
packets_buffer = buffer(packets, m)';               % Group bits into bits per symbol
sym_idx = bi2de(packets_buffer, 'left-msb')'+1;     % Bits to symbol index
x = s(sym_idx);                                     % Look up symbols using the indices  
x_upsample = upsample(x, fsfd);                     % Space the symbols fsfd apart, to enable pulse shaping using conv.
x_upsample(end-fsfd+2:end) = [];

end

