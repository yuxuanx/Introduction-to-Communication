function [ signal_modulated ] = baseband2passband( pulse_tr_RC_samp, f_carrier, fs )
% Move the tranmitted signal to the carrier
t=0:1/fs:(length(pulse_tr_RC_samp)-1)/fs;
% Modulation
Icarrier = sqrt(2)*(real(pulse_tr_RC_samp)).*cos(2*pi*f_carrier*t);     
Qcarrier = sqrt(2)*(imag(pulse_tr_RC_samp)).*sin(2*pi*f_carrier*t);
signal_modulated=Icarrier+Qcarrier;
signal_modulated = signal_modulated./abs(max(signal_modulated)); % Normalization          
% N = length(signal_modulated);
% P = fftshift(fft(signal_modulated,N));                       % Fourier transform
% fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
% figure;
% plot(fvec,20*log10(abs(P)));

end

