% Test the findPoints script

% Test data
n_points = 1000;
t = 1:n_points;
phase = (0.1) * sin(2 * pi * t / n_points);
shifts = [110, 230, 410, 650, 850];
n_shifts = length(shifts);

for i=1:n_shifts
    phase(shifts(i):end) = phase(shifts(i):end) + 1;
end

points_per_bin = 20;

points = find_points(phase, points_per_bin);

% There are 5 shifts and they should occur in bins 6, 12, 21, 33, 43

find(points == 1)