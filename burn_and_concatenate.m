function [points, history] = burn_and_concatenate(initial_points, break_index, bins_per_window, n_windows)
% Create the history matrix for the point process 
% 
% Args:
%   points (array): Point process of interest
%   bins_per_window (int): Size of each history window
%   n_windows (int): Number of history windows for each channel
%
% Returns:
%   history (array, n_bins by n_parameters): Point process history

n_bins = size(initial_points, 2);
n_channels = size(initial_points, 1);

n_segments = length(break_index);
break_index = [break_index, n_bins];

n_parameters = 1 + n_channels * n_windows;

burn_length = n_windows * bins_per_window;

bins_per_segment = zeros(1, n_segments);
for segment = 1:n_segments
    bins_per_segment(segment) = max(break_index(segment + 1) - break_index(segment) - burn_length, 0);
end
total_bins = sum(bins_per_segment);

points = zeros(n_channels, total_bins);
history = zeros(total_bins, n_parameters);

pointer = 1;
for segment = 1:n_segments 
    if bins_per_segment(segment) > 0
        segment_history = make_history(initial_points(:, break_index(segment):(break_index(segment + 1) - 1)), bins_per_window, n_windows);
        segment_history(1:burn_length, :) = [];
        segment_points = initial_points(:, (break_index(segment) + burn_length):(break_index(segment + 1) - 1));
        points(:, pointer : pointer + bins_per_segment - 1) = segment_points;
        history(pointer:pointer + bins_per_segment - 1, :) = segment_history;
    end
end

