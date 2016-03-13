function [points]=findPoints(Phase, delta, alpha)
% Find the phase shift events for each channel in a dataset. Time series is
% binned into equal size bins of size delta.
% Args:
%   Phase (array): Instant phase time series for each signal.
%   delta (int): Temporal resolution. Number of points in each bin.
%   alpha (float): Significance level for shift identification
% 
% Returns:
%   dN (array): Bins with phase shift events for each signal.

% Number of channels
n_channels = size(Phase, 1);

% Number of bins
n_bins = floor(size(Phase, 2) / delta);

% Number of events in each bin
points = zeros(n_channels, n_bins);

%Partition the time series into bins of length delta
interval = (1:n_bins) * delta;

% Identify points for each channel
for c = 1:n_channels   
    % Identify latency of phase shift events using phase derivative    
    shift_times = shift_identification_PD(Phase(c, :), alpha);
    % Assign shifts to dN
    points(c, 1) = sum(shift_times < interval(1));
    for i=2:n_bins
        points(c, i) = sum(shift_times < interval(i)) - points(c, i-1);
    end
end