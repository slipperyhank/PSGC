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
n_bins = 40000;
bins_per_window = 50;
n_windows = 2;
p_event = 0.02;
points_per_bin = 4;
frequency_band = [8, 10];
frequency = 9;
n_points = n_bins * points_per_bin;
max_windows = 3;
shift_magnitude = 2;
alpha = 0.01;
signal2noise_ratio = 0.5;
to_burn = 0;
break_index = 1;

data = zeros(n_channels, n_points);

t = 1:n_points;
sampling_rate = 250;

% Events in channel 1 excite channel 2 in the 6-10 bin lag 
channel1_parameters = [log(p_event), -4, 0, 0, 0];
channel2_parameters = [log(p_event), 0, 2, -6, 0];
model_parameters = {channel1_parameters, channel2_parameters};

n_iterations = 100;

result_psgc = zeros(n_channels);
result_order = zeros(1, n_channels);

n_failed = 0;

for iteration = 1:n_iterations
    iteration
    [points, history] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window, to_burn);
    for channel = 1:n_channels
        shift_index = find(points(channel, :) == 1) * points_per_bin;
        phase = zeros(1, n_points);
        phase(shift_index) = shift_magnitude;
        phase = cumsum(phase);
        data(channel, :) = sin(2 * pi * frequency * t / sampling_rate + phase);
        data(channel, :) = data(channel, :) / norm(data(channel, :));
        noise = randn(size(data(channel, :)));
        noise = noise / norm(noise);
        data(channel, :) = signal2noise_ratio * data(channel, :) + (1 - signal2noise_ratio) * noise;
    end
    try
        [pval_matrix, parameter_estimates, model_order] = psgc_analysis(data, break_index, frequency_band, sampling_rate, points_per_bin, bins_per_window, max_windows, alpha);
    catch
        n_failed = n_failed + 1;
        pval_matrix = zeros(n_channels);
        model_order = zeros(1, n_channels);
    end
    result_psgc = result_psgc + (pval_matrix < 0.05);
    result_order = result_order + model_order;
end

result_psgc = result_psgc / (n_iterations - n_failed);
result_order = result_order / (n_iterations - n_failed);


% In 100 iterations: Power 100% to identify the three true connections,
% false positive rate of 4% on the connection that does not exist.

