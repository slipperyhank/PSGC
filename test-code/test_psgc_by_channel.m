% Randomly simulate point process from a specific model. Test the 
% psgc_by_channel function to test the LRTs power and false positive rate.

% Default values:
% min_step_size = 0.01
% n_channels = 2
% n_points = 10000
% bins_per_window = 5
% n_window = 2
% p_event = 0.005

n_channels = 2;
n_bins = 10000;
bins_per_window = 5;
n_windows = 2;
p_event = 0.005;


% Events in channel 1 excite channel 2 in the 6-10 bin lag 

channel1_parameters = [log(p_event), 0, 0, 0, 0];
channel2_parameters = [log(p_event), 0, 2, 0, 0];
model_parameters = {channel1_parameters, channel2_parameters};

n_iterations = 100;
pvalues = zeros(1, n_channels);

for iteration = 1:n_iterations
    iteration
    [points, history] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window);
    try
        [pvalues(iteration, :)] = psgc_by_channel(points(2, :), history, n_windows);
    catch
        pvalues(iteration, :) = -1;
    end
end

% In a 100 randomly simulated datasets, there was a 4% false positive rate
% and power of 96% respectively



n_channels = 2;
n_bins = 40000;
bins_per_window = 5;
n_windows = 2;
p_event = 0.01;
max_windows = 3;

channel1_parameters = [log(p_event), 0, 0, 0, 0];
channel2_parameters = [log(p_event), -1, 0, 0, 0];
model_parameters = {channel1_parameters, channel2_parameters};

n_iterations = 100;
pvalues = zeros(1, n_channels);

for iteration = 1:n_iterations
    iteration
    [points, history] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window);
    try
        [pvalues(iteration, :)] = psgc_by_channel(points(2, :), history, n_windows);
    catch
        pvalues(iteration, :) = -1;
    end
end

% With n_bins = 10000, p_event = 0.01 and 100 iterations, power is 31% and false positive
% rate is 4%

% With n_bins = 25000, p_event = 0.01 and 100 iterations, power is 56% and
% false postive rate is 5%

% With n_bins = 40000, p_event = 0.01 and 100 iterations, power is 89% and
% false positive rate is 4% 

% There is reduced power to find inhibitory singals - probably because
% p_event is already so close to 0. 


n_channels = 2;
n_bins = 10000;
bins_per_window = 5;
n_windows = 2;
p_event = 0.005;

channel1_parameters = [log(p_event), 0, 0, 0, 1.5];
channel2_parameters = [log(p_event), 0, 1.5, 0, 0];
model_parameters = {channel1_parameters, channel2_parameters};

n_iterations = 100;
pvalues = zeros(1, n_channels);

for iteration = 1:n_iterations
    iteration
    [points, ~] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window);
    history = make_history(points, bins_per_window, n_windows + 2);
    try
        [pvalues(iteration, :)] = psgc_by_channel(points(2, :), history, n_windows + 2);
    catch
        pvalues(iteration, :) = -1;
    end
end

% With n_bins = 10000, and 100 iterations, power is 76% and false positive
% rate is 4%

% With n_bins = 10000, and 100 iterations and wrong window size (4 instead
% of 2, makes it harder to decouple), power is 66% and false positive rate
% is 9%.

% With n_bins = 25000, and 100 iterations, power is 99% and false positive
% rate is 8%