% Randomly simulate point process from a specific model. Reconstruct time
% series, add noise, and test find_points. 

n_channels = 2;
n_bins = 5000;
bins_per_window = 50;
n_windows = 2;
p_event = 0.02;
points_per_bin = 4;
frequency_band = [8, 10];
frequency = 9;
n_points = n_bins * points_per_bin;
max_windows = 3;
shift_magnitude = 2;
band_centre = frequency;
band_width = frequency_band(2) - band_centre;
data = zeros(n_channels, n_points);
alpha = 0.01;
t = 1:n_points;
sampling_rate = 250;
signal2noise_ratio = 0.5;
error_threshold = 5;

n_iterations = 100;
true_positive_rate = zeros(1, n_iterations);
false_positive_rate = zeros(1, n_iterations);

% Events in channel 1 excite channel 2 in the 6-10 bin lag 
channel1_parameters = [log(p_event), -4, 0, 0, 0];
channel2_parameters = [log(p_event), 0, 2, -6, 0];  
model_parameters = {channel1_parameters, channel2_parameters};
phase = zeros(n_channels, n_points);

for iteration = 1:n_iterations
    [point_process] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window);

    for channel = 1:n_channels
        shift_index = find(point_process(channel, :) == 1) * points_per_bin;
        true_phase = zeros(1, n_points);
        true_phase(shift_index) = shift_magnitude;
        true_phase = cumsum(true_phase);
        data(channel, :) = sin(2 * pi * frequency * t / sampling_rate + true_phase);
        data(channel, :) = data(channel, :) / norm(data(channel, :));
        noise = randn(size(data(channel, :)));
        noise = noise / norm(noise);
        data(channel, :) = signal2noise_ratio * data(channel, :) + (1 - signal2noise_ratio) * noise;
        phase(channel, :) = instant_phase(data(channel, :), sampling_rate, band_centre, band_width);
    end

    points = find_points(phase, points_per_bin);

    true_points_channel1 = find(point_process(1, :) > 0);
    estimated_points_channel1 = find(points(1, :) > 0);
    n_points_channel1 = length(true_points_channel1);
    accuracy_channel1 = zeros(1, n_points_channel1);
    estimated_point_index1 = zeros(1, n_points_channel1);
    for i = 1:n_points_channel1
        [accuracy_channel1(i), estimated_point_index1(i)] = min(abs(true_points_channel1(i) - estimated_points_channel1));
    end
    estimated_point_index1(accuracy_channel1 >= error_threshold) = '';

    true_points_channel2 = find(point_process(2, :) > 0);
    estimated_points_channel2 = find(points(2, :) > 0);
    n_points_channel2 = length(true_points_channel2);
    accuracy_channel2 = zeros(1, n_points_channel2);
    estimated_point_index2 = zeros(1, n_points_channel2);
    for i = 1:n_points_channel2
        [accuracy_channel2(i), estimated_point_index2(i)] = min(abs(true_points_channel2(i) - estimated_points_channel2));
    end
    estimated_point_index2(accuracy_channel2 >= error_threshold) = '';

    true_positives = length(unique(estimated_point_index1)) + length(unique(estimated_point_index2));
    false_positives = length(estimated_points_channel1) + length(estimated_points_channel2) - true_positives;

    true_positive_rate(iteration) = true_positives / (n_points_channel1 + n_points_channel2);
    false_positive_rate(iteration) = false_positives / (length(estimated_points_channel1) + length(estimated_points_channel2));
end

mean(true_positive_rate)
mean(false_positive_rate)

% Default values 
% bins_per_window = 50;
% n_windows = 2;
% p_event = 0.02;
% points_per_bin = 4;
% frequency_band = [8, 10];
% frequency = 9;
% max_windows = 3;
% shift_magnitude = 2;
% alpha = 0.01;
% sampling_rate = 250;
% signal2noise_ratio = 0.5
% channel1_parameters = [log(p_event), -4, 0, 0, 0];
% channel2_parameters = [log(p_event), 0, 2, -6, 0];
% error_threshold = 5

%%%%%%%%%%%
% Results %
%%%%%%%%%%%

% In the default settings over 5 trials
% TPR = 0.9516
% FPR = 0.0318

% Decrease alpha to 0.001
% TPR = 0.9536
% FPR = 0.0250

% The result does not seem highly dependent on alpha at this scale. Suspect
% with lower shift magnitude or greater SNR, alpha becomes more relevant.

% FPR decreases with alpha, as expected.

% Change channel parameters so there are more closely spaced events
% channel2_parameters = [log(p_event), 0, 2, -2, 0];

% In this case, the algorithm for identifying true positives can be flawed:
%  a single estimated event can be marked as satisfying two true events

% TPR = 0.7223
% FPR = 0.0838

% why does the FPR go up? maybe more true shifts are just mislocalized?
% same channel parameters with increased error threshold = 7

% TPR = 0.7583
% FPR = 0.0381




% Other things to check
% power decreases with shift magnitude
% power decreases (slowly) with SNR
% power decreases when more shift occur close together
% 
