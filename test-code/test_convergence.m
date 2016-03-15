% Randomly generate phase signals with specific overall probably of event
% Generate history, test for PSGC convergence. 

% Run with min_step_size = 0.01
% Default parameters are
% n_channels = 3
% n_points = 1000
% bins_per_window = 5
% n_window = 2

n_channels = 3;
n_points = 1000;
bins_per_window = 2;
n_window = 2;

n_converge = zeros(1, 20);

for i=1:20
    p_event = i / 40;
    points = rand(n_channels, n_points);
    points = points < p_event;
    history = make_history(points, bins_per_window, n_window);
    likelihood = zeros(1, n_channels);
    for channel = 1:n_channels
        try
            [likelihood(channel), final_parameters] = fit_model(points(channel, :), history);
        catch
            likelihood(channel) = 0;
        end
    end
    n_converge(i) = sum(likelihood ~= 0);
end


% With the default settings, p <= 0.25 is the approximate criteria
% for convergence.

% With 2 bins_per_window, p <= 0.5 is the approximate criteria
% for convergence.

% With 4 n_windows, p <= 0.1 is the approximate criteria

% With n_channels = 5, p <= (0.15, 0.175) is the approximate criteria

% With n_channels = 5 and n_window = 4, p <= 0.075 is the approximate
% criteria

% Problems can crop up quickly if history windows have more than one event
% in them