function P = instant_phase(x, sRate, centre, width)
% Estimate the instantaneous phase of a time series in frequency band using
% the Hilbert transform.
%
% Args: 
%   x (array): time series to be analyzed
%   sRate (int): Sampling rate of the signal
%   centre (float): Centre of the frequency band of interest
%   width (float): Width of the band of interest (must be symmetric)
% 
% Returns:
%   P (array): Estimated instantaneous phase of the signal.

% Coerce x to a double
x = double(x);

% Transform frequency component to zero
t = 0:(length(x)-1);
Y1 = x .* sin(2 * pi * centre / sRate * t);
Y2 = x .* cos(2 * pi * centre / sRate * t);

% Low-pass filter power outside band of interest
[b,a]=butter(4, 2 * width / sRate, 'low');
F1 = filtfilt(b, a, Y1);        
F2 = filtfilt(b, a, Y2);    

% Solve for phase
P = atan2(F2, F1);

