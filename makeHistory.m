function history = makeHistory(points, bins_per_window, n_window)
% Create the history matrix for the point process 
% 
% Args:
%   points (array): Point process of interest
%   bins_per_window (int): Size of each history window
%   n_windows (int): Number of history windows for each channel
%
% Returns:
%   history (array, n_bins by n_parameters): Point process history

n_bins = size(points, 2);
n_channels = size(points, 1);

% Initialize history
history = zeros(n_bins, (1 + n_channels * n_window));

% Intercept term
history(:, 1) = 1;

% First n_window * bins_per_window bins are initial conditions
for bin = (n_window * bins_per_window + 1):n_bins
    for window = 1:n_window
        for channel = 1:n_channels
            % In the history, windows for each channel are grouped together
            history(bin, 1 + (channel - 1) * n_window + window) = sum(points(channel, (bin - window * bins_per_window):(bin - 1 - (window - 1) * bins_per_window)));       
        end
    end
end
          
    
end