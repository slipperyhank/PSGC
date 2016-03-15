function [model_order, best_parameters, AIC] = findModel(points, bins_per_window, max_windows)
% Find the model order for each channel
% Args:
%   points (array): The number of shifts in each bin
%   bins_per_window (int): The number of data points in each history window
%   max_windows (int): The maximum number of history windows to try
%
% Returns:
    % model_order (int): The optimal number of history windows
    % best_parameters (array): The set of parameters for the optimal model
    

% Initialize model order, maximum likelihood and parameters values
n_channels = size(points, 1);
model_order = zeros(1, n_channels);
AIC = zeros(1, n_channels);
best_parameters = cell(1, n_channels);

for n_windows = 1:max_windows;
    n_windows
    % Calculate H for each channel
    history = makeHistory(points, bins_per_window, n_windows);
    % Calculate likelihood and parameter estimates for each
    % channel
    for channel = 1:n_channels
        channel
        [likelihood, parameters] = fit_model(points(channel, (n_windows * bins_per_window +1):end), history((n_windows * bins_per_window + 1):end, :));
        % If this likelihood is greater than all previous, update the AIC,
        % model_order and best_parameters
        if 2 * (1 + n_channels * n_windows - likelihood) < AIC(channel)
            AIC(channel) = 2 * (1 + n_channels * n_windows - likelihood);
            best_parameters{channel} = parameters;
            model_order(channel) = n_windows;
        end
    end
end