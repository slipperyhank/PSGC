function [points, history] = burn_and_concatenate(initial_points, bin_boundary_markers, bins_per_window, n_windows)
% Create the history matrix for the point process 
% 
% Args:
%   points (array): Point process of interest
%   bin_boundary_markers (array(int)): Array of bin indices, marking the
%      start of a new segment. Final value should be n_bins + 1.
%   bins_per_window (int): Size of each history window
%   n_windows (int): Number of history windows for each channel
%
% Returns:
%   points (array, n_channels by n_bins): New point process with initial
%      values burned away. 
%   history (array, n_bins by n_parameters): Point process history

% If bin_bounardy_markers does not mark the end of the file, add it
if bin_boundary_markers(end) ~= size(initial_points, 2) + 1
    bin_boundary_markers = [bin_boundary_markers, size(initial_points, 2) + 1];
end

% Initialize from inputs
n_channels = size(initial_points, 1);
n_segments = length(bin_boundary_markers) - 1;
n_parameters = 1 + n_channels * n_windows;
burn_length = n_windows * bins_per_window;

bins_per_segment = zeros(1, n_segments);
for segment = 1:n_segments
    bins_per_segment(segment) = bin_boundary_markers(segment + 1) - bin_boundary_markers(segment) - burn_length;
end
total_bins = sum(bins_per_segment);

points = zeros(n_channels, total_bins);
history = zeros(total_bins, n_parameters);

pointer = 1;
for segment = 1:n_segments 
    if bins_per_segment(segment) > 0
        segment_history = make_history(initial_points(:, bin_boundary_markers(segment):(bin_boundary_markers(segment + 1) - 1)), bins_per_window, n_windows);
        segment_history(1:burn_length, :) = [];
        segment_points = initial_points(:, (bin_boundary_markers(segment) + burn_length):(bin_boundary_markers(segment + 1) - 1));
        points(:, pointer : pointer + bins_per_segment(segment) - 1) = segment_points;
        history(pointer:pointer + bins_per_segment(segment) - 1, :) = segment_history;
        pointer = pointer + bins_per_segment(segment);
    end
end

