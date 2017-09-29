% Test the burn_and_concatenate function

% Test data
n_bins = 100;
n_channels = 2;
n_windows = 2;
bins_per_window = 5;
bin_boundary_markers = [1, 51, 101];
n_segments = length(bin_boundary_markers) - 1;

point_process = zeros(2, n_bins);
point_process(1, [9, 21, 31, 66, 76, 86]) = 1;
point_process(2, [4, 26, 46, 56, 61, 71]) = 1;

% There are 6 events in channel 1
% The events occur in original bins [9, 21, 31, 66, 76, 86]
% burn_and_concatenate removes bins [1:10, 51:60] and after reindex the
% remaining 5 events should be in [11, 21, 46, 56, 66]

% There are 6 events in channel 2
% The events occur in original bins [4, 26, 46, 56, 61, 71]
% burn_and_concatenate removes bins [1:10, 51:60] and after reindex the
% remaining 4 events should be in [16, 36, 41, 51]
 
% The size of the history should be 80 bins (10 burned from each segment): 
% Column 1 is all ones
% Column 2 has original events at (10-14, 22-26, 32-36, 67-71, 77-81, 87-91)
% then reindex after burn removal (1-4, 12-16, 22-26, 47-51, 57-61, 67-71)
% Column 3 has original events (15-19, 27-31, 37-41, 72-76, 82-86, 92-96)
% then reindex after burn removal (5-9, 17-21, 27-31, 52-56, 62-66, 72-76)

% Column 4 has original events at (5-9, 27-31, 47-51, 57-61, 62-66, 72-76)
% then reindex after burn removal (17-21, 37-40, 41, 42-46, 52-56)
% Column 5 has original events at (10-14, 32-36, 52-56, 62-66, 67-71, 77-81)
% then reindex after burn removal (1-4, 22-26, 42-46, 47-51, 57-61)

% For channel two they should occur in bins 4, 13, 15

[points, history] = burn_and_concatenate(point_process, bin_boundary_markers, bins_per_window, n_windows);

isequal(find(points(1, :) == 1), [11, 21, 46, 56, 66])
isequal(find(points(2, :) == 1), [16, 36, 41, 51])

isequal(sum(history(:, 1)), n_bins - n_breaks * n_windows * bins_per_window)
isequal(find(history(:, 2) == 1), [1:4, 12:16, 22:26, 47:51, 57:61, 67:71]')
isequal(find(history(:, 3) == 1), [5:9, 17:21, 27:31, 52:56, 62:66, 72:76]')
isequal(find(history(:, 4) == 1), [17:21, 37:46, 52:56]')
isequal(find(history(:, 5) == 1), [1:4, 22:26, 42:51, 57:61]')

