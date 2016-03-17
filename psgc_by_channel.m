function [pvalues, statistics, reduced_likelihood] = psgc_by_channel(points, history, n_windows, varargin)
% Find the PSGC statistics and p-values for incoming connections to a
% specific channel
%
% Args:
%   points (array(int)): The point process for the channel of interest
%   history (array(int)): The CIF history at the appropriate model order
%   n_windows (int): The model order for the channel of interest
%   varargin:
%   {1} parameters (array(float)): Option estimated parameters that have 
%                                   been pre-computed for the full model   
%
% Returns:
%   pvalues (array(float)): pvalues testing for significant PSGC for each 
%                           input channel
%   statistics (array(float)): Chi^2 values for likelihood ratio test
%   reduced_likelihood (array(float)): Likelihood value of reduced models

n_bins = length(points);
n_parameters = size(history, 2);
n_channels = (n_parameters - 1) / n_windows;

% Output variables
reduced_likelihood = zeros(1, n_channels);

% Optional initial conditions
if nargin < 4
    [likelihood, parameters] = fit_model(points, history); 
elseif nargin == 4
    parameters = varargin{1};
    if ~iscolumn(parameters)
        parameters = parameters';
    end
    [likelihood, ~] = fit_model(points, history, parameters); 
else
    error('Must be at least 2 inputs. Option input is estiamted parameters of full model')
end

for channel = 1:n_channels
    % Create reduced history with the contribution of channel removed
    % TODO - test reduce code
    reduced_history = history;
    reduced_history(:, ((channel - 1) * n_windows + 2):(channel * n_windows + 1)) = [];
    reduced_parameters = parameters;
    reduced_parameters(((channel - 1) * n_windows + 2):(channel * n_windows + 1)) = [];
    % fit the reduced point process model
    [reduced_likelihood(channel), ~] = fit_model(points, reduced_history, reduced_parameters);                 
end

statistics = -2 * (reduced_likelihood - likelihood);
pvalues = 1 - chi2cdf(statistics, n_windows);



