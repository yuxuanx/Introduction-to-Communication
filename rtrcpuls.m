function [y] = rtrcpuls(a,tau,fs,span)
% This function is a modified version of the one given in the text, Program 2.2-2
% a: Roll off/ Excess bandwidth factor
% tau: Symbol time 
% fs: sampling frequency at which the continuous pulse is sampled
% span: Defines width for truncation. Number of tau periods on either side of the peak of the pulse
% The pulse has a one sided bandwidth, BW = (1+alpha)/(2*tau);

t_positive = eps:(1/fs):span*tau;  % Replace 0 with eps (smallest +ve number MATLAB can produce) to prevent NANs
t = [-fliplr(t_positive(2:end)) t_positive];
tpi = pi/tau; amtpi = tpi*(1-a); aptpi = tpi*(1 + a);
ac = 4*a/tau; at = 16*a^2/tau^2;
y = (sin(amtpi*t) + (ac*t).*cos(aptpi*t))./(tpi*t.*(1-at*t.^2));
norm_factor = sqrt(sum(y.^2));
y = y/norm_factor; % Normlaize the pulse to have unit energy