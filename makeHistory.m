function H = makeHistory(dN, W, M)
% Create the history matrix for the point process dN
% 
% Args:
%   dN (array): Point process of interest
%   W (int): Size of each history window
%   M (int): Number of history windows
%
% Returns:
%   H (array): Matrix of point process history

K = size(dN, 2);
nChan = size(dN, 1);

% Initialize history
H = zeros(K, (1 + nChan * M));

% Intercept term
H(:, 1) = 1;

for i = (M * W + 1):K
    for m = 1:M
        for c = 1:nChan
            H(i, 1 + c * (M-1) + m) = sum(dN(c, (i - m * W):(i - 1 - (m - 1) * W)));       
        end
    end
end
          
    
end