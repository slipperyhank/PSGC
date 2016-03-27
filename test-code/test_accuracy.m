% Randomly simulate point process from a specific model. Test the fit_model
% function to see if parameter estimates converge to the correct values. 

% Default values:
% min_step_size = 0.01
% n_channels = 2
% n_points = 1000
% bins_per_window = 5
% n_window = 2
% p_event = 0.005

n_channels = 2;
n_bins = 50000;
bins_per_window = 5;
n_windows = 2;
p_event = 0.005;
to_burn = 1;

% First two independent processes
channel1_parameters = [log(p_event), 0, 0, 0, 0];
channel2_parameters = [log(p_event), 0, 0, 0, 0];
model_parameters = {channel1_parameters, channel2_parameters};
n_parameters = length(channel1_parameters);

n_iterations = 50;
parameters = cell(n_iterations, n_channels);
likelihood = zeros(1, n_channels);
estimated_parameters = cell(1, n_channels);
for channel = 1:n_channels
    estimated_parameters{channel} = zeros(1, n_parameters);
end

% Simulate process and history
for iteration = 1:n_iterations
    [points, history] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window, to_burn);
    for channel = 1:n_channels
        [likelihood(channel), parameters{iteration, channel}] = fit_model(points(channel, :), history);
    end
    for channel = 1:n_channels
        parameters{iteration, channel}(parameters{channel} == -Inf) = 0;
        estimated_parameters{channel} = estimated_parameters{channel} + parameters{iteration, channel} / n_iterations;
    end
end
    
% Rate is estimate with good accuracy. Other parameters appear to be highly
% influenced by random chance - high variance. Need to test significance 
% to determine if this is a problem 




% Events in channel 1 increase first excite, then inhibit events in channel
% 2. Channel 2 has no effect on channel 1. 
channel1_parameters = [log(p_event), 0, 0, 0, 0];
channel2_parameters = [log(p_event), 2, -1, 0, 0];
model_parameters = {channel1_parameters, channel2_parameters};
n_parameters = length(channel1_parameters);
n_iterations = 5;
parameters = cell(n_iterations, n_channels);
likelihood = zeros(1, n_channels);
estimated_parameters = cell(1, n_channels);
for channel = 1:n_channels
    estimated_parameters{channel} = zeros(1, n_parameters);
end

% Simulate process and history
for iteration = 1:n_iterations
    [points, history] = simulate_from_model(model_parameters, n_bins, n_windows, bins_per_window);
    for channel = 1:n_channels
        try
            [likelihood(channel), parameters{iteration, channel}] = fit_model(points(channel, :), history);
        catch
            likelihood(channel) = 0;
            parameters{iteration, channel} = zeros(size(channel1_parameters));
        end
    end
    for channel = 1:n_channels
        parameters{iteration, channel}(parameters{iteration, channel} == -Inf) = 0;
        estimated_parameters{channel} = estimated_parameters{channel} + parameters{iteration, channel} / n_iterations;
    end
end

% Non-zero parameters estimated with reasonable accuracy. Zero parameters
% appear to be highly variable. Need to test for significance.


% Crude statistical analysis of estimated parameters

nboot = 5000;
bootfun = @(x)mean(x);

for channel = 1:n_channels
    for parameter = 1:n_parameters
        y = zeros(1, n_iterations);
        for iteration = 1:n_iterations
            y(iteration) = parameters{iteration, channel}(parameter);
        end
        ci = bootci(nboot, {bootfun, y}, 'type', 'per');
        [model_parameters{channel}(parameter), ci(1), ci(2)]
    end
end

% At n_bins = 10000 and p_event = 0.002, the ci's for the 0 parameters
% were all positive.

% At n_bins = 50000 and p_event = 0.005, the true parameters falls in the
% ci for many of the parameters, and is very close in the rest.

% Most parameters appear towards the border of their CI. Suspect there is
% some bias in the results, due to finding inhibitory and excitatory
% signals in the data when none exist in the model. Should happen less with
% more data points. 

% This bias appear to have dimenished after baseline bins burned away
% TODO: recheck the above more thoroughly 

