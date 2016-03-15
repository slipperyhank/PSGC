% Test the findPoints script

% Test data
n_points = 100;
t = 1:n_points;
phase1 = (0.1) * sin(2 * pi * t / n_points);
phase2 = phase1;

shifts1 = [12, 54, 66];
shifts2 = [16, 42, 72];

n_shifts = length(shifts1);

for i=1:n_shifts
    phase1(shifts1(i):end) = phase1(shifts1(i):end) + 1;
    phase2(shifts2(i):end) = phase2(shifts2(i):end) + 1;
end

phase = [phase1; phase2];

alpha = 0.05;
bin_size = 5;
n_bins = floor(n_points / bin_size);

points = find_points(phase, bin_size, alpha);

% There are 3 for each channel 
% For channel one they should occur in bins 3, 11, 14
% For channel two they should occur in bins 4, 9, 15

find(points(1, :) == 1)
find(points(2, :) == 1)

bins_per_window = 1;
n_windows = 2;

history = make_history(points, bins_per_window, n_windows);

% The first column is an intercept
% In the second column, there should be ones at (4, 12, 15)
% In the third column, there should be ones at (5, 13, 16)
% In the fourth column, there should be ones at (5, 10, 16)
% In the fifth column, there should be ones at (6, 11, 17)
 
isequal(sum(history(:, 1)), n_bins)
isequal(find(history(:, 2) == 1), [4; 12; 15])
isequal(find(history(:, 3) == 1), [5; 13; 16])
isequal(find(history(:, 4) == 1), [5; 10; 16])
isequal(find(history(:, 5) == 1), [6; 11; 17])

