function [pval_matrix, parameter_estimates, statistic_matrix, model_order] = psgc_analysis(data, break_index, frequency_band, sampling_rate, points_per_bin, bins_per_window, possible_models, alpha)
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
'Calculating instant phase'
for channel = 1:n_channels
    data(channel, :) = instant_phase(data(channel, :), sampling_rate, band_centre, band_width);
end

% Identify phase shift events
'Transforming into point process'
point_process = zeros(n_channels, n_bins);
for channel = 1:n_channels
    point_process(channel, :) = find_points(data(channel, :), points_per_bin, alpha);
end

% Identify optimal model order for each channel
'Finding model orders'
parameter_estimates = cell(1, n_channels);
flag = zeros(1, n_channels);
if length(possible_models) == 1
    model_order = zeros(1, n_channels) + possible_models;
    [points, history] = burn_and_concatenate(point_process, break_index, bins_per_window, possible_models);
    for channel = 1:n_channels
        try
            [~, parameter_estimates{channel}] = fit_model(points(channel, :), history);
        catch
            flag(channel) = 1;
        end
    end
elseif length(possible_models) == n_channels
    model_order = possible_models;
    unique_orders = unique(model_order);
    n_models = length(unique_orders);
    for model = 1:n_models
        [points, history] = burn_and_concatenate(point_process, break_index, bins_per_window, unique_orders(model));        
        for channel = 1:n_channels
            try
                [~, parameter_estimates{channel}] = fit_model(points(channel, :), history);
            catch
                flag(channel) = 1;
            end 
        end
    end
elseif length(possible_models) == 2
    [model_order, parameter_estimates] = find_model(point_process, break_index, bins_per_window, possible_models(1):possible_models(2));
end

% Perform likelihood ratio tests
'Calculating PSGC'
failed_channels = [];
pval_matrix = zeros(n_channels);
statistic_matrix = zeros(n_channels);
for channel = 1:n_channels
    channel
    if model_order(channel) > 0 && flag(channel) == 0
        try
            [points, history] = burn_and_concatenate(point_process, break_index, bins_per_window, model_order(channel));
            [pval_matrix(channel, :), statistic_matrix(channel, :)] = psgc_by_channel(points(channel, :), history, model_order(channel), parameter_estimates{channel});
        catch
            failed_channels = [failed_channels, channel];
        end
    elseif iscell(model_order)
        pval_matrix(channel, :) = 1;
        statistic_matrix(channel, :) = 0;
    end
end






