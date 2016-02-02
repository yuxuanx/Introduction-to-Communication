  function transmitter(packet,fc)
% packet: information bits to be tranmistted
% fc: carrier frequency
% Parameters below like "rb","rs" can be caculated by the fomular BW =
% (1+alpha)/(2*tau) as well. 
% In the situation below the BW = (1+0.4)/(2*1/rs) = 210;
fs = 12e3;                                   % sampling frequency [Hz]
rb = 600;                                    % bit rate [bit/sec]
M = 4;                                       % Number of symbols in the constellation (QPSK, M=4)
m = log2(M);                                 % Number of bits per symbol
rs = rb/m;                                   % Symbol rate
fsrs = fs/rs;                                % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity)
packet = packet';
barkerBits = [0 0 0 0 0 1 1 0 0 1 0 1 0];
markerBits = repmat(barkerBits, 1, 16);      % Duplicate barkerbits to 16 times
load('syncBits.mat')                         % 68 bits synchonization bits
dataBits = [markerBits, syncBits, packet];

RC_puls = rtrcpuls(0.4,1/rs,fs,6);           % Pulse shaping function
x_upsample = bits2symbols(dataBits, fsrs, m);   % Bits to symbols
pulse_tr_RC_samp = pulseShaping(RC_puls, x_upsample, fsrs, fs); % Pulse shaping
signal_modulated = baseband2passband(pulse_tr_RC_samp ,fc, fs); % Modulation
sound(signal_modulated,fs);                  % Play the transmitted signal
audiowrite('trial.wav',signal_modulated,fs)
end