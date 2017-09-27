function [boundary_markers] = eeglab_get_boundary_markers(EEG, varargin)
% Find all boundary markers from and eeglab EEG file
% Input: 
%    EEG: an eeglab EEG data structure
%    event_name (string; optional): name of boundary events 
%                                   default is 'boundary'
% Output: 
%    array (int): Index of each bounary markers, signifying the beginning
%                 of a new epoch

if nargin == 1
    event_name = 'boundary';
else
    event_name = varargin{1};
end

n_events = length(EEG.event);
n_boundaries = 0;
boundary_markers = 1;

for i = 1:n_events
    if strcmp(EEG.event(i).type, event_name)
        n_boundaries = n_boundaries + 1;
        boundary_markers(n_boundaries) = int64(EEG.event(i).latency + 0.5);
    end
end

