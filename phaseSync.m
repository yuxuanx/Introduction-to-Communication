function [ mf_phase ] = phaseSync( syncSymbol, mf_downsample )
sumArg = 0;
conjSync = conj(syncSymbol);
for k =1:length(conjSync)
    arg = angle(mf_downsample(k)*conjSync(k));
    sumArg = sumArg + arg;
end
phihat = sumArg/length(conjSync);
mf_phase = mf_downsample * exp(-1j*phihat);


end

