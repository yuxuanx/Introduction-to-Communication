function [ mf_samp ] = matchedFilter( RC_puls, Icarrier_remove, Qcarrier_remove, fsfd)
% Make the transmitted baseband signal go through a matched filter for
% decision making
MF_puls=fliplr(RC_puls);                            % Matched filter is a time-reversed pulse shaping filter 
mf=conv(MF_puls,Icarrier_remove+1j*Qcarrier_remove);% Make the signal through the matched filter
mf_samp = mf(fsfd*6:end-fsfd*5);
% eyed.fsfd=fsfd;
% eyed.r=mf_samp;
% eyediagram(eyed.r, eyed.fsfd);                      % Plot the eye diagram
% figure;
% plot(real(mf_samp));
% N = length(mf_samp);
% P_mf = fftshift(fft(mf_samp,N));                    % Fourier transform
% fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
% figure;
% plot(fvec,20*log10(abs(P_mf)));
end

