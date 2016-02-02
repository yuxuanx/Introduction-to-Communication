function [ Xhat ] = symbols2bits( Ifinal, Qfinal, mf_phase )
% Transfer the symbols back to data bits
final=[Ifinal(1:end)',Qfinal(1:end)'];
finalbits=zeros(length(mf_phase),2);
for i=1:length(final)
    if final(i,1)==1 && final(i,2)==1
        finalbits(i,1)=0;
        finalbits(i,2)=0;
    elseif final(i,1)==1 && final(i,2)==-1
        finalbits(i,1)=0;
        finalbits(i,2)=1;
    elseif final(i,1)==-1 && final(i,2)==-1
        finalbits(i,1)=1;
        finalbits(i,2)=0;
    else
        finalbits(i,1)=1;
        finalbits(i,2)=1;
    end
end
finalbits=finalbits';
Xhat=reshape(finalbits,1,length(final)*2);
end

