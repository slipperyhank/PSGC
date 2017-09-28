function [pval_matrix, parameter_estimates, statistic_matrix, model_order] = psgc_analysis(data, boundary_markers, sampling_rate, frequency_band, points_per_bin, bins_per_window, maximum_model_order, varargin)
% Perform phase shift granger causality on set of signals.
%
% Args: 
%   data (n_channels by n_points): A matrix of data. 
%   boundary_markers (array(int)): Array of boundary indices. For
%      contiguous data, boundary_markers = 1.
%   sampling_rate (int): Number of samples per second.
%   frequency_band (array): Frequency band of interest [lower, upper].
%   points_per_bin (int): Number of data points in each observtion bin.
%   bins_per_window (int): Number of observation bins per history window.
%   maximum_model_order (int): Maximum number of history windows to be
%      considered in model selection. 
%   alpha (float): Threshold for identifying phase shift events.
%   reverse (boolean): True if data should be time reversed.
%
% Returns:
%   pval_matrix (n_channels by n_channels): Matrix of p-values testing the
%      hypothesis that row i does not predict column j.
%   parameter_estimates (array): The parameters values for the model of 
%       each channels point process. 
%   statistic_matrix ():
%   model_order ():


if nargin == 8
    reverse = varargin{1};
else
    reverse = false;
end

% Number of channels
n_channels = size(data, 1);

% Find frequency band centre and width
band_centre = (frequency_band(2) + frequency_band(1)) / 2;
band_width = frequency_band(2) - band_centre;

% Add marker to the end of the data
boundary_markers = [boundary_markers, size(data, 2) + 1];

% Minimum segment length needed to contribute to PSGC model
min_points = maximum_model_order * bins_per_window * points_per_bin;

% Number of segments
n_segments = sum(diff(boundary_markers) > min_points);

% Set the end points, and number of bins, for each segment 
segment_endpoints = zeros(n_segments, 2);
bins_per_segment = zeros(1, n_segments);
count = 0;
for i = 1:(length(boundary_markers) - 1)
    % Determine size of potential segment
    n_points = boundary_markers(i + 1) - boundary_markers(i);
    if n_points > min_points
        count = count + 1;
        segment_endpoints(count, 1) = boundary_markers(i);
        segment_endpoints(count, 2) = boundary_markers(i + 1) - 1;
        bins_per_segment(count) = floor(n_points / points_per_bin);
    end
end

% Set bin boundary markers
bin_boundary_markers = cumsum(bins_per_segment) + 1;
bin_boundary_markers = [1, bin_boundary_markers];

% Total number of bins
n_total_bins = sum(bins_per_segment);

% Initialize point process
point_process = zeros(n_channels, n_total_bins);

for segment = 1:n_segments   
    % Grab specific segment data for all channels
    segment_data = data(:, segment_endpoints(segment, 1):segment_endpoints(segment, 2));
    % Flip data if time-reversed anlaysis is desired
    if reverse
        data = fliplr(data);
    end
    % Calculate the instantaneous phase for each signal
    for channel = 1:n_channels
        segment_data(channel, :) = instant_phase(segment_data(channel, :), sampling_rate, band_centre, band_width);
    end
    % Identify phase shift events
    for channel = 1:n_channels
        point_process(channel, bin_boundary_markers(segment):(bin_boundary_markers(segment + 1) - 1)) = find_points(segment_data(channel, :), points_per_bin);
    end
end


% Identify optimal model order for each channel
'Finding model orders'
[model_order, parameter_estimates] = find_model(point_process, bin_boundary_markers, bins_per_window, maximum_model_order);

% Perform likelihood ratio tests
'Calculating PSGC'
failed_channels = [];
pval_matrix = zeros(n_channels);
statistic_matrix = zeros(n_channels);
for channel = 1:n_channels
    channel
    try
        [points, history] = burn_and_concatenate(point_process, bin_boundary_markers, bins_per_window, model_order(channel));
        [pval_matrix(channel, :), statistic_matrix(channel, :)] = psgc_by_channel(points(channel, :), history, model_order(channel), parameter_estimates{channel});
    catch
        failed_channels = [failed_channels, channel];
        pval_matrix(channel, :) = 1;
        statistic_matrix(channel, :) = 0;
    end
end







