function [ Icarrier_remove, Qcarrier_remove ] = passband2baseband( signal_modulated, f_carrier, fs )
% Move the modulated signal to baseband to demodulate it
t=0:1/fs:(length(signal_modulated)-1)/fs;
Icarrier_remove=sqrt(2)*signal_modulated.*cos(2*pi*f_carrier*t);
Qcarrier_remove=sqrt(2)*signal_modulated.*sin(2*pi*f_carrier*t);
% signal_demodulated=Icarrier_remove+Qcarrier_remove;
% N = length(signal_demodulated);
% P = fftshift(fft(signal_demodulated,N));                % Fourier tranform
% fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
% figure;
% plot(fvec,20*log10(abs(P)));
end

