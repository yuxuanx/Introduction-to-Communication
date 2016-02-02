function [ errorrate ] = errorComputing( packets, Xhat )
% Compute the error rate 
diff=packets-Xhat;
error=find(diff~=0);
errorrate=length(error)/length(final);

end

