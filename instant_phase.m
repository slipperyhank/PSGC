function [phase, magnitude] = instant_phase(signal, sampling_rate, band_centre, band_width)
% Estimate the instantaneous phase of a time series in frequency band using
% the Hilbert transform.
%
% Args: 
%   signal (array): time series to be analyzed
%   sampling_rate (int): Sampling rate of the signal
%   band_centre (float): Centre of the frequency band of interest
%   band_width (float): Width of the band of interest (must be symmetric)
% 
% Returns:
%   phase (array): Estimated instantaneous phase of the signal.
%   magnitude (array): Estimated instantaneous phase of the signal.

% TODO - make filter order an optional input
filter_order = 4;

% Coerce x to a double
signal = double(signal);

% Transform frequency component to zero
t = 0:(length(signal) - 1);

in_phase = signal .* sin(2 * pi * band_centre / sampling_rate * t);
quadrature = signal .* cos(2 * pi * band_centre / sampling_rate * t);

% Low-pass filter power outside band of interest
[b, a] = butter(filter_order, 2 * band_width / sampling_rate, 'low');

% Filter the in_phase and in_quadrature terms
in_phase = filtfilt(b, a, in_phase);        
quadrature = filtfilt(b, a, quadrature);    

% Solve for phase
phase = atan2(quadrature, in_phase);
magnitude = sqrt(in_phase .^ 2 + quadrature .^ 2);

