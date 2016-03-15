% Randomly simulate point process from a specific model. Test the find_model
% function to see if model orders are correct. 

% Default values:
% min_step_size = 0.01
% n_channels = 2
% n_points = 10000
% bins_per_window = 5
% n_window = 2
% p_event = 0.005
% max_windows = 5

n_channels = 2;
n_bins = 40000;
bins_per_window = 5;
n_windows = 2;
p_event = 0.005;
max_windows = 3;

% First two independent processes
channel1_parameters = [log(p_event), 0, 0, 0, 0];
channel2_parameters = [log(p_event), 0, 0, 0, 0];
model_parameters = {channel1_parameters, channel2_parameters};

n_iterations = 5;
parameters = cell(1, n_iterations);
model_order = zeros(n_iterations, n_channels);
AIC = zeros(n_iterations, n_channels);

% Simulate process and history
for iteration = 1:n_iterations
    iteration
    [points, history] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window);
    [model_order(iteration, :), parameters{iteration}, AIC(iteration, :)] = find_model(points, bins_per_window, max_windows);
end
    
% Model order estimate is good - all correct in 5 iterations.
% In more simulations it occasionally errors by 1, and rarely by 2.


% Events in channel 1 increase first excite, then inhibit events in channel
% 2. Channel 2 has no effect on channel 1. 
channel1_parameters = [log(p_event), 0, 0, 0, 0];
channel2_parameters = [log(p_event), 2, -1, 0, 0];
model_parameters = {channel1_parameters, channel2_parameters};

n_iterations = 5;
parameters = cell(1, n_iterations);
model_order = zeros(n_iterations, n_channels);
AIC = zeros(n_iterations, n_channels);

% Simulate process and history
for iteration = 1:n_iterations
    iteration
    [points, history] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window);
    [model_order(iteration, :), parameters{iteration}, AIC(iteration, :)] = find_model(points, bins_per_window, max_windows);
end

% With n_bins = 10000, and n_iterations = 10, average error is 0.5 and
% 12/20 correct, 6/20 off by one, 2/20 off by 2, 0/20 off by three
% (possible because max_windows = 3).

% With n_bins = 40000, and n_iterations = 5, average error is 0.3 and
% 7/10 correct, 3/10 off by one, 0/10 off by two or three.
