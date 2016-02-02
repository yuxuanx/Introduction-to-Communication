function [ mf_sync ] = symbolSync( syncSymbol, fsrs, mf_samp )
sum = zeros(1,10000);
for tsamp=1:6000
    for k=1:length(syncSymbol)
        sum(tsamp)=sum(tsamp)+mf_samp((k-1)*fsrs+tsamp)*conj(syncSymbol(k));
    end
end
[~,tsamp]=max(abs(sum));
mf_sync = mf_samp(tsamp:end);


end

