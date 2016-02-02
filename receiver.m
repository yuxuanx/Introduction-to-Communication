function [Xhat, psd, const, eyed] = receiver(tout,fc)

fs = 12e3;                                   % sampling frequency [Hz]
rb = 600;                                    % bit rate [bit/sec]
M = 4;                                       % Number of symbols in the constellation (QPSK, M=4)
m = log2(M);                                 % Number of bits per symbol
rs = rb/m;                                   % Symbol rate
fsrs = fs/rs;                                % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity)

load('syncBits.mat')
load('syncSymbol.mat')
load('marker_modulated.mat')
load('RC_puls.mat')
tic
signal_modulated = signalRecording(rs, fs ,marker_modulated);
timeElapsed = toc;
if timeElapsed < tout
    [Icarrier_remove, Qcarrier_remove] = passband2baseband(signal_modulated, fc, fs);
    mf_samp = matchedFilter( RC_puls, Icarrier_remove, Qcarrier_remove, fsrs);
    mf_sync = symbolSync( syncSymbol, fsrs, mf_samp );
    mf_downsample = downsample(mf_sync, fsrs);          % Downsampling the signal after matched filter
    mf_phase = phaseSync( syncSymbol, mf_downsample );
    [ Ifinal, Qfinal ] = decisionMaking( mf_phase );
    Xhat = symbols2bits( Ifinal, Qfinal, mf_phase );
    bitsRestore = frameSync( Xhat, syncBits );
    info_samp = mf_phase(length(syncSymbol)+1:length(syncSymbol)+216);
    Xhat = bitsRestore(length(syncBits)+1:end); % received information bits
    [pvalue,fvalue] = pwelch(signal_modulated,window,[],2048,fs); % received signal
    pvalue = pvalue/max(pvalue);
    pvalue = 10*log10(pvalue);
    psd = struct('p',pvalue,'f',fvalue);
    const = info_samp;  % after downsampling
    eyed = [info_samp, fsrs];
else
    Xhat = [];
    psd = struct('p',0,'f',0);
    const = [];
    eyed = [];
end
end