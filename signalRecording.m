function [ signal_modulated ] = signalRecording(rs, fs ,marker_modulated)
cor=0;
threshold = 10;  
recObj = audiorecorder(fs,8,1);              % Set record object

while max(cor)<threshold
   recordblocking(recObj,7/rs);
   signal_modulated = getaudiodata(recObj);
   cor=xcorr(signal_modulated, marker_modulated);
end
recordblocking(recObj,350/rs);
signal_modulated = getaudiodata(recObj);
signal_modulated=signal_modulated';
end

