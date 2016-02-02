%% bits to Symbols
clc; clear all; close all

% headBits = [0 1 0 1 0 1 0 1];
% trainingLen = 60;
% trainingBits = randsrc(1,trainingLen,[0 1]);
fs = 12e3;                                          % sampling frequency [Hz]
rb = 720;                                           % bit rate [bit/sec]
infoLen = 432;

% Len = length(headBits) + trainingLen + infoLen;
% number of bits to transmit
% Constellation or bit to symbol mapping
s = [exp(0) exp(1j*pi/4) exp(1j*3*pi/4) exp(1j*pi/2) exp(1j*7*pi/4) exp(1j*3*pi/2) exp(1j*pi) exp(1j*5*pi/4)]; % Constellation 1 - QPSK/4-QAM
                                                    % s = exp(1i*((0:3)*pi/2 + pi/4));
M = length(s);                                      % Number of symbols in the constellation
m = log2(M);                                        % Number of bits per symbol
fd = rb/m;                                          % Symbol rate
fsfd = fs/fd;                                       % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity)
rng('shuffle')
infoBits = randsrc(1,infoLen,[0 1]);
% Information bits
barkerBits = [0 0 0 0 0 1 1 0 0 1 0 1 0];
markerBits = repmat(barkerBits, 1, 18);
syncBits = randsrc(1,90,[0,1]);
Bits = [syncBits, infoBits];
dataBits = [markerBits, syncBits, infoBits];
b_buffer = buffer(dataBits, m)';                           % Group bits into bits per symbol
sym_idx = bi2de(b_buffer, 'left-msb')'+1;           % Bits to symbol index
x = s(sym_idx);                                     % Look up symbols using the indices  
x_upsample = upsample(x, fsfd);                     % Space the symbols fsfd apart, to enable pulse shaping using conv.
x_upsample(end-fsfd+2:end) = [];
%% RRC pulse
span = 6;
beta = 0.4;
RC_puls = rtrcpuls(beta,1/fd,fs,span);         % Root raised cosine pulse shaping filter
pulse_tr_RC = conv(RC_puls,x_upsample);
pulse_tr_RC_samp = pulse_tr_RC(fsfd*6:end-fsfd*5);  % Discard the last fsfd*(span-1) data
N = length(pulse_tr_RC_samp);
P = fftshift(fft(pulse_tr_RC_samp,N));              % Fourier transform
fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
figure; plot(fvec,20*log10(abs(P)));
%% Carrier
f_carrier=3000;                                     % Carrier frequency
t=0:1/fs:(length(pulse_tr_RC_samp)-1)/fs;
% Modulation
Icarrier = sqrt(2)*(real(pulse_tr_RC_samp)).*cos(2*pi*f_carrier*t);     
Qcarrier = sqrt(2)*(imag(pulse_tr_RC_samp)).*sin(2*pi*f_carrier*t);
carrier=Icarrier+Qcarrier;
N = length(carrier);
P = fftshift(fft(carrier,N));                       % Fourier transform
fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
figure;
plot(fvec,20*log10(abs(P)));
carrier = carrier./abs(max(carrier));
% marker_frequency = 8e3;
modulateBits = [barkerBits barkerBits];
marker_upsample = bits2symbols(barkerBits, fsfd, m);
pulse_marker = pulseShaping(RC_puls, marker_upsample, fsfd, fs);
marker_modulated = baseband2passband(pulse_marker ,f_carrier, fs);
marker_modulated = marker_modulated./abs(max(marker_modulated));
% sound(marker_modulated,fs);                                  % Play the transmitted signal
sound(carrier,fs);  
% audio = [marker_modulated, carrier];
% sound(audio,fs);
audiowrite('trial.wav',carrier,fs)
%% Signal through AWGN
% snr=100;                                             % Signal-to-noise ratio
% carrier_noise=awgn(Icarrier+1j*Qcarrier,snr);       % Through Gussain white noise channel
% N = length(carrier_noise);
% fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
% P_noise = fftshift(fft(carrier_noise,N));           % Fourier transform
% figure;
% plot(fvec,20*log10(abs(P_noise)));
cor=0;
threshold = 10;
recObj = audiorecorder(fs,8,1);              % Set record object
% recordblocking(recObj,tout);
% audioData = getaudiodata(recObj);
while max(cor)<threshold
   record(recObj,9/fd);
   stop(recObj);
   signal_modulated = getaudiodata(recObj);
%    signal_modulated = signalRecording(10/fd, fs);
   cor=xcorr(signal_modulated, marker_modulated);
end
figure();
plot(abs(cor));
recordblocking(recObj,360/fd);
carrier_noise = getaudiodata(recObj);
figure();
plot(carrier_noise);
% envelope = abs(hilbert(carrier_noise));
% figure;
% plot(envelope); % used for detect where voice signal begins
%% Remove carrier
% Demodulation
% Icarrier_remove=sqrt(2)*real(carrier_noise).*cos(2*pi*f_carrier*t);
% Qcarrier_remove=sqrt(2)*imag(carrier_noise).*sin(2*pi*f_carrier*t);
t=0:1/fs:((length(carrier_noise)-1)/fs);
Icarrier_remove=sqrt(2)*carrier_noise'.*cos(2*pi*f_carrier*t);
Qcarrier_remove=sqrt(2)*carrier_noise'.*sin(2*pi*f_carrier*t);
carrier_remove=Icarrier_remove+Qcarrier_remove;
N = length(carrier_remove);
P = fftshift(fft(carrier_remove,N));                % Fourier tranform
fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
figure;
plot(fvec,20*log10(abs(P)));
%% Matched filter
MF_puls=fliplr(RC_puls);                            % Matched filter is a time-reversed pulse shaping filter 
figure;
plot(MF_puls);
mf=conv(MF_puls,Icarrier_remove+1j*Qcarrier_remove);% Make the signal through the matched filter
mf_samp = mf(fsfd*6:end-fsfd*5);
eyed.fsfd=fsfd;
eyed.r=mf_samp;
eyediagram(eyed.r, eyed.fsfd);                      % Plot the eye diagram
N = length(mf_samp);
P_mf = fftshift(fft(mf_samp,N));                    % Fourier transform
fvec = (fs/N)*(-floor(N/2):1:ceil(N/2)-1);
figure;
plot(fvec,20*log10(abs(P_mf)));
figure;
plot(real(mf_samp));
%% Symbol Synchronization
packets_buffer = buffer(syncBits, m)';               % Group bits into bits per symbol
idxx = bi2de(packets_buffer, 'left-msb')'+1;     % Bits to symbol index
syncSymbol = s(idxx); 
sum = zeros(1,10000);
for tsamp=1:6000
    for k=1:length(syncSymbol)
        sum(tsamp)=sum(tsamp)+mf_samp((k-1)*fsfd+tsamp)*conj(syncSymbol(k));
    end
end
[~,tsamp]=max(abs(sum));
mf_samp2 = mf_samp(tsamp:end);
figure;
plot(real(mf_samp2));
eyed.fsfd=fsfd;
eyed.r=mf_samp2;
eyediagram(eyed.r, eyed.fsfd); 
% [acor,lag] = xcorr(mf_samp,syncSymbol);
% [~,I] = max(abs(acor));
% timeDiff = lag(I);
%% Decision making (Sample at Ts)
mf_downsample = downsample(mf_samp2, fsfd);          % Downsampling the signal after matched filter
mf_downsample = mf_downsample(1:end-1);
% mf_downsample = mf_samp(timeDiff:fsfd:timeDiff+fsfd*(length(dataBits)/2-1));
scatterplot(mf_downsample); grid on;                % Plot the constellation of the signal after downsampling

[psd.p,psd.f]=pwelch(mf_downsample);                
figure;
plot(psd.f,20*log10(psd.p));                        % Plot the power spectral density

% Phase Synchronization
sumArg = 0;
conjSync = conj(syncSymbol);
for k =1:length(conjSync)
    arg = angle(mf_downsample(k)*conjSync(k));
    sumArg = sumArg + arg;
end
phihat = sumArg/length(conjSync);
mf_downsample = mf_downsample * exp(-1j*phihat);
scatterplot(mf_downsample); grid on; 
% threshold=0;                                      % Decide the threshold of to make decision
% realpart=real(mf_downsample);
% imagpart=imag(mf_downsample);
% % Decision making progress
% for i=1:length(mf_downsample)
%     if realpart(i)>=threshold
%         Ifinal(i)=1;
%     else
%         Ifinal(i)=-1;
%     end
%     if imagpart(i)>=threshold
%         Qfinal(i)=1;
%     else
%         Qfinal(i)=-1;
%     end
% end
Ifinal=zeros(1,length(mf_downsample));
Qfinal=zeros(1,length(mf_downsample));
    for i=1:length(mf_downsample)
        D1=norm(s(1)-mf_downsample(i));%Calculate euclidean distance to each
        D2=norm(s(2)-mf_downsample(i));%point of our constellation
        D3=norm(s(3)-mf_downsample(i));
        D4=norm(s(4)-mf_downsample(i));
        D5=norm(s(5)-mf_downsample(i));
        D6=norm(s(6)-mf_downsample(i));
        D7=norm(s(7)-mf_downsample(i));
        D8=norm(s(8)-mf_downsample(i));
        D=[D1 D2 D3 D4 D5 D6 D7 D8]; %Put all the distance in one vector
        [~, I]=min(D);   %Search for the index of the smallest value
        Ifinal(i)=real(s(I)); %And use this index to determine which symbol was send
        Qfinal(i)=imag(s(I));
    end
%% Symbols to bits
final=Ifinal(1:end)+1j*Qfinal(1:end);
finalbits=zeros(length(mf_downsample),3);
for i=1:length(final)
    if final(i)==s(1)
        finalbits(i,1)=0;
        finalbits(i,2)=0;
        finalbits(i,3)=0;
    elseif final(i)==s(2)
        finalbits(i,1)=0;
        finalbits(i,2)=0;
        finalbits(i,3)=1;
    elseif final(i)==s(3)
        finalbits(i,1)=0;
        finalbits(i,2)=1;
        finalbits(i,3)=0;
    elseif final(i)==s(4)
        finalbits(i,1)=0;
        finalbits(i,2)=1;
        finalbits(i,3)=1;
    elseif final(i)==s(5)
        finalbits(i,1)=1;
        finalbits(i,2)=0;
        finalbits(i,3)=0;
    elseif final(i)==s(6)
        finalbits(i,1)=1;
        finalbits(i,2)=0;
        finalbits(i,3)=1;
    elseif final(i)==s(7)
        finalbits(i,1)=1;
        finalbits(i,2)=1;
        finalbits(i,3)=0;
    else
        finalbits(i,1)=1;
        finalbits(i,2)=1;
        finalbits(i,3)=1;
    end
end
finalbits=finalbits';
Xhatt=reshape(finalbits,1,length(final)*3);
%% Frame Synchonization
% corr = conv(Xhatt,fliplr(syncBits));
% [tmp, idx] = max(corr);
% Xhat = Xhatt(1+idx:length(Bits)+idx);

c=zeros(1,length(syncBits)+1);
for m=0:length(syncBits)-1
    for i=1:length(syncBits)
    c(m+1) = c(m+1) + syncBits(i)*Xhatt(i+m);
    end
end
[~,idx]=max(c);
Xhat = Xhatt(idx:idx+length(Bits)-1);
%% Calculating the error rate
diff=Bits-Xhat;
error=find(diff~=0);
errorrate=length(error)/length(Xhat);
