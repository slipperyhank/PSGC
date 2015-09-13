function [N1,c] = shift_identification_PD(P, alpha)
% Identify the phase shift events in a time series
% of phase values.
%
% Args:
%   P (array) - time series of phase values
%   alpha (float) - significance level
%
% Returns:
%   N1 (array) - point process of phase shifts
%   c (float) - actual PD threshold used


% Unwrap the phase if not already done
P = unwrap(P);

% Number of time points
N = length(P);

% Estimate phase derivative and subtract mean
dP = (P(3:end) - P(1:(end - 2))) / 2;
dP = dP - mean(dP);

% Critical values is the 1 - alpha / 2 quantile of the Normal distribution
crit = norminv(1 - alpha / 2);
sigma = std(dP);
c = sigma * crit;

% Continue removing shift values and updating the threshold until
% convergence
while 1==1
    ind = abs(dP) > c;
    sigma = min(std(dP(ind == 0)), sigma);
    if c == sigma * crit;
        break
    end
    c = sigma * crit;
end

ind = [0, ind, 0];
dP = [0, dP, 0];

lower = find(ind(1:(end - 1)) - ind(2:end) == -1) + 1;
upper = find(ind(1:(end - 1)) - ind(2:end) == 1);

nShifts = length(lower);

N1 = zeros(1, nShifts);
for i=1:nShifts
    N1(i) = find(abs(dP) == max(abs(dP(lower(i):upper(i)))), 1);
end









        