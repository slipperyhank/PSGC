function [dN]=findPoints(Phase, delta, alpha)
% Find the phase shift events for each channel in a dataset.
% Args:
%   Phase (array): Instant phase time series for each signal.
%   delta (int): Temporal resolution. Number of points in each bin.
%   alpha (float): Significance level for shift identification
% 
% Returns:
%   dN (array): Bins with phase shift events for each signal.

nChan = size(Phase, 1);

N = size(Phase, 2);

dN = zeros(nChan, N);

%Partition the time series into K bins of length delta
K=floor(N / delta);
interval = 1 + (1:K) * delta;

for c = 1:nChan    
    % Identify phase shift events using PD    
    N1 = shift_identification_PD(Phase(c, :), alpha);
    % Assign shifts to dN
    dN(1) = sum(N1 < interval(1));
    for i=2:K
        dN(i) = sum(N1 < interval(i)) - dN(i-1);
    end
end
        

