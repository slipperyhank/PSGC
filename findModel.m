function [maxM, g, AIC] = findModel(dN, W, nM)
% Find the model order for each channel
% Args:
%   dN (array): The number of shifts in each bin
%   W (int): The number of data points in each history window
%   nM (int): The maximum number of history windows to try
%
% Returns:
    % maxM (int): The optimal number of history windows
    % g (array): The set of parameters for the optimal model
    

% Initialize model order, maximum likelihood and parameters values
nChan = size(dN, 1);
maxM = zeros(1, nChan);
AIC = zeros(1, nChan);
g = cell(1, nChan);

for M=1:nM;
     % Calculate H for each channel
     H = makeHistory(dN, M, W, delta);
     % Calculate likelihood and parameter estimates for each
     % channel
     for c = 1:nChan
         [Likelihood, Gamma]=PSGC_noIC(dN(c, :), H);
         % If this likelihood is better than all previous, update the AIC,
         % maxM and g
         if 2 * (1 + nChan * M - Likelihood) < AIC(c)
             AIC(c) = 2 * (1 + nChan * M - Likelihood);
             g{c} = Gamma;
             maxM(c) = M;
         end
     end
end