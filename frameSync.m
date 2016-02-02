function [ pack ] = frameSync( Xhat, syncBits )
c=zeros(1,length(syncBits)+1);
for m=0:length(syncBits)-1
    for i=1:length(syncBits)
    c(m+1) = c(m+1) + syncBits(i)*Xhat(i+m);
    end
end
[~,idx]=max(c);
pack = Xhat(idx:(idx+length(syncBits)+432-1));

end

