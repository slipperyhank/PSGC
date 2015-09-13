function P = instant_phase(x, sRate, freq, width)

% Estimate the instantaneous phase of a time series x
% In the band centred at freq with length 2*width

% Phase is not unwrapped, and no Burn period is remove

x = double(x);
 
t = 0:(length(x)-1);

[b,a]=butter(4, 2 * width / sRate, 'low');

Y1 = x .* sin(2 * pi * freq / sRate * t);
Y2 = x .* cos(2 * pi * freq / sRate * t);

F1 = filtfilt(b, a, Y1);        
F2 = filtfilt(b, a, Y2);    

P = atan2(F2, F1);

