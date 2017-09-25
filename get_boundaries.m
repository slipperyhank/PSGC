function [boundaries] = get_boundaries(EEG)
% TODO add option input string for name of boundary event
% For now the default is only 'boundary'
% Input: an eeglab EEG data structure
% Output: an array of bounary indices. Each index signifies the beginning
% of a new epoch

n_events = length(EEG.event);
n_boundaries = 0;
boundaries = 0;

for i = 1:n_events
    if strcmp(EEG.event(i).type, 'boundary')
        n_boundaries = n_boundaries + 1;
        boundaries(n_boundaries) = EEG.event(i).latency + 0.5;
    end
end

if boundaries(1) == 0
    boundaries = 1;
end
