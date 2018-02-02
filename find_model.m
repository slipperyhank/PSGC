function [best_models, best_parameters, AIC] = find_model(point_process, bin_boundary_markers, bins_per_window, maximum_model_order)
% Find the model order for each channel
% Args:
%   point_process (array: n_channel by n_bins): The number of shifts in 
%      each bin for each channel.
%   bin_bounary_markers (array(int)): Indices of the bins that mark the
%      start of a new segment. 
%   bins_per_window (int): The number of data points in each history
%      window.
%   max_windows (int): The maximum number of history windows to try during
%      model selection. 
%
% Returns:
%   model_order (int): The optimal number of history windows
%   best_parameters (array): The set of parameters for the optimal model
%   AIC (array(float)):
    
% Initialize model order, maximum likelihood and parameters values
n_channels = size(point_process, 1);
best_models = zeros(1, n_channels);
AIC = zeros(1, n_channels) + Inf;
best_parameters = cell(1, n_channels);

for model_order = 0:maximum_model_order;
    model_order
    % Calculate history for this value of n_windows
    [points, history] = burn_and_concatenate(point_process, bin_boundary_markers, bins_per_window, model_order);
    n_parameters = size(history, 2);
    % Calculate likelihood and parameter estimates for each
    % channel
    for channel = 1:n_channels
        try
            [likelihood, parameters] = fit_model(points(channel, :), history);
        catch
            likelihood = -Inf;
            parameters = zeros(1, n_parameters);
        end
        % If this likelihood is greater than all previous, update the AIC,
        % model_order and best_parameters
        if 2 * (1 + n_channels * model_order - likelihood) < AIC(channel)
            AIC(channel) = 2 * (1 + n_channels * model_order - likelihood);
            best_parameters{channel} = parameters;
            best_models(channel) = model_order;
        end
    end
end