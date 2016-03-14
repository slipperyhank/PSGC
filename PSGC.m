function [likelihood, parameters]=PSGC(points, history, bin_time, varargin)
% Numerical solve for the maximum likelihood estimates of a PSGC model 
% using Newtons method
%
% Args: 
%   points (array): The binned phase shift events a single channel
%   history (array): The history values of each channel
%   bin_time (float): The time in seconds of each bin
%   initial (array(float)): Option intial conditions for optimization
%
% Returns:
%   L (float): Maximum likelihood value
%   Gamma (array(float)): Maximum likelihood estimates

% Input points should be a column vector. If it is a row, transpose
if isrow(points)
    points = points';
end
% Validate that points is now a column
if ~iscolumn(points)
    error('points input must be a column vector')
end

% First remove all parameters that are trivially -Inf
% i.e., if this event has occured, then the current process will
% not have an event
inhibitory_parameters = find(points'*history == 0);
[inhibitory_bins, ~] = find(history(:, inhibitory_parameters) == 1);

history(:, inhibitory_parameters) = [];
history(inhibitory_bins, :) = [];
points(inhibitory_bins) = [];

% Now remove all parameters taht are trivially +Inf
% i.e., if this event has occured, the the current process will
% definitely have an event.

excitatory_parameters = find(points'*history == sum(history));
[excitatory_bins, ~] = find(history(:, excitatory_parameters) == 1);

history(:, excitatory_parameters) = [];
history(excitatory_bins, :) = [];
points(excitatory_bins) = [];

% Flag for determining recursion or model fit error or no problems
flag=0;

% convergence criteria
epsilon=0.00001;

% number of segments
n_bins = length(points);

%number of parameters
n_parameters = size(history, 2);

%Initialize output
parameters = zeros(n_parameters, 1);
likelihood = 0;

% Optional initial condition input
% Initialize parameters (X)
if nargin < 4
    X = zeros(n_parameters, 1);
elseif nargin == 4
    X = varargin{1};
else
    error('Wrong number of inputs, must be three or four')
end

% Index for parameters
index = 1:n_parameters;

% Bad measures parameters which go below the minimum (-10)
%Bad = X;

%initialize Jacobian (J) and Derivative (F) and intensity (lambda) for
%optimzation method
J = zeros(n_parameters);
F = zeros(n_parameters, 1);
lambda = zeros(n_bins, 1);

% Optimize 
while 1 == 1
    for bin = 1:n_bins
        % Calculate lambda
        lambda(bin) = exp(sum(X' .* history(bin, :)));
        % If probability of a shift grows too large, model fails
        if lambda(bin) * bin_time > 1
            flag=1
        end
    end
    lambda * bin_time
    if flag==1
        error('Probability greater than 1')
        break
    end
    
    % Calculate gradiant (F) and Jacobian (J) for newtons method
    for parameter1 = 1:n_parameters
        F(parameter1) = -sum(points .* history(:, parameter1) + history(:, parameter1) .* (bin_time * lambda .* ...
            (points - 1)) ./ (1 - bin_time .* lambda));
        for parameter2 = 1:n_parameters
            J(parameter1, parameter2) = sum(history(:, parameter1) .* history(:, parameter2) .* (points - 1) .* ...
                bin_time .* lambda ./ (1 - lambda .* bin_time) .^ 2);
        end
    end
    % Calculate change in parameters dX. Magnitude is reduced to 20% to
    % prevent parameters from over shooting and jumping outside the support 
    % of the model
    F
    J
    dX = J \ F
    X = X + dX / 200
    % If any of the parameters are below -10, go to recursion step
    if min(X) < -100
        %Bad = X < -10;
        flag = 2
        break;
    end
    % when convergence is reached, exit loop
    
    if max(abs(dX)) < epsilon 
        break
    end
end

% flag=0: No errors, calculate output
if flag == 0
    parameters = X;
    for bin = 1:n_bins
        lambda(bin) = exp(sum(X' .* history(bin,:))); 
    end
    likelihood = sum((points .* log(lambda .* bin_time) + (1 - points) .* log(1 - lambda .* bin_time)));
end

% flag=1: error, must be looked at
if flag == 1
    %likelihood = 0;
    %parameters = X .* 0;
    flag
end

% flag=2: parameters below -10, remove from optimization and called
% function again
if flag == 2
    flag
    %parameters(index(Bad==1)) = -10;
    %index(Bad==1) = '';
    %ROld = '';
    %[likelihood, temp] = PSGC2(points, R, Bad, ROld, X(Bad==0));
    %parameters(index) = temp;
end
