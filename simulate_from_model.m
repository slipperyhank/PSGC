function [points, history] = simulate_from_model(parameter_cell, n_bins, n_windows, bins_per_window)
% Simulate from the binomial approximation to a Poisson model with CIF 
%
% Args: 
%   n_windows (int): The number of history windows to be used
%   bins_per_window (int): Number of bins per history window
%   parameters (cell(array(float))): The parameters of the fit model. Each
%   cell entry is the parameter set for a specific channel

% Returns:
%   points (array(int)): The simulated point process for each channel
%   history (array(int)): The corresponding history values


% need to set M and gamma for the specific Model
% Use max M in the head, and then zero out other parameters
n_parameters = length(parameter_cell{1});
n_channels = (n_parameters - 1) / n_windows;

if ~rem(n_channels, 1) == 0 
    error('Wrong parameter structure')
end

history = zeros(n_bins, n_parameters);
history(:, 1) = 1;
points = zeros(n_channels, n_bins);

% Set initial conditions without history
for channel = 1:n_channels
    for bin = 1:(n_windows * bins_per_window)
        if rand() < exp(parameter_cell{channel}(1))
            points(channel, bin) = 1;
        end
    end
end

for bin = (n_windows * bins_per_window + 1):n_bins
   for channel = 1:n_channels
       for window = 1:n_windows
           history(bin, 1 + (channel - 1) * n_windows + window) = sum(points(channel, (bin - window * bins_per_window):(bin - 1 - (window - 1) * bins_per_window)));
       end
   end
   for channel = 1:n_channels
       parameters = parameter_cell{channel};
       if rand() < exp(sum(parameters.*history(bin, :)))
           points(channel, bin) = 1;
       end
   end
end
  