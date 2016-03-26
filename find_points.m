function [points] = find_points(phase, points_per_bin, alpha)
% Find the phase shift events for each channel in a dataset. Time series is
% binned into equal size bins of size delta.
% Args:
%   Phase (array): Instant phase time series for each signal.
%   bin_width (int): Temporal resolution. Number of points in each bin.
%   alpha (float): Significance level for shift identification
% 
% Returns:
%   points (array): Point process - list of bins with phase shift events 
%                   for each channel.


% Number of channels
n_channels = size(phase, 1);

% Number of bins
n_bins = floor(size(phase, 2) / points_per_bin);

% Number of events in each bin
points = zeros(n_channels, n_bins);

%Partition the time series into bins of length delta
interval = (1:n_bins) * points_per_bin;

% Identify points for each channel
for channel = 1:n_channels    
    % Identify latency of phase shift events using phase derivative    
    shift_times = shift_identification(phase(channel, :), alpha);
    % Assign shifts to dN
    points(channel, 1) = sum(shift_times < interval(1));
    for i=2:n_bins
        points(channel, i) = sum(shift_times < interval(i)) - sum(shift_times < interval(i-1));
    end
end