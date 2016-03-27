function [pval_matrix, parameter_estimates, model_order] = psgc_analysis(data, break_index, frequency_band, sampling_rate, points_per_bin, bins_per_window, max_windows, alpha)
% Perform phase shift granger causality on set of signals.
%
% Args: 
%   data (n_channels by n_points): A matrix of data. Default assumption is contiguous.  
%   frequency_band (array): Frequency band of interest [lower, upper]
%   bins_per_window (int): Number of bins per history window
%
% Returns:
%   pval_matrix (n_channels by n_channels): Matrix of p-values testing the
%      hypothesis that row i does not predict column j.
%   parameter_estimates (array): The parameters values for the model of 
%       each channels point process. 
 
n_channels = size(data, 1);
n_points = size(data, 2);
n_bins = floor(n_points / points_per_bin);



% Find frequency band centre and width
band_centre = (frequency_band(2) + frequency_band(1)) / 2;
band_width = frequency_band(2) - band_centre;

% Calculate the instantaneous phase for each signal
for channel = 1:n_channels
    data(channel, :) = instant_phase(data(channel, :), sampling_rate, band_centre, band_width);
end

% Identify phase shift events
point_process = zeros(n_channels, n_bins);
for channel = 1:n_channels
    point_process(channel, :) = find_points(data(channel, :), points_per_bin, alpha);
end

% Identify optimal model order for each channel
[model_order, parameter_estimates] = find_model(point_process, break_index, bins_per_window, max_windows);

% Perform likelihood ratio tests
pval_matrix = zeros(n_channels);
for channel = 1:n_channels
    if model_order(channel) > 0
        [points, history] = burn_and_concatenate(point_process, break_index, bins_per_window, model_order(channel));
        pval_matrix(channel, :) = psgc_by_channel(points(channel, :), history, model_order(channel), parameter_estimates{channel});
    else
        pval_matrix(channel, :) = 1;
    end
end




