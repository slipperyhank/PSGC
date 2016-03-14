function [likelihood, final_parameters]=PSGC(points, history, varargin)
% Numerical solve for the maximum likelihood estimates of a PSGC model 
% using Newtons method
%
% Args: 
%   points (array): The binned phase shift events a single channel
%   history (array): The history values of each channel
%   bin_time (float): The time in seconds of each bin
%   varargin:
%   {1} initial (array(float)): Option intial conditions for optimization
%   {2} step_size (float): Step size for newton optimization
%
% Returns:
%   L (float): Maximum likelihood value
%   Gamma (array(float)): Maximum likelihood estimates
%   indices (cell(array(int))): List of inhibitory and excitatory
%                               parameters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process varargin and validate inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create initial versions incase need to recall function
initial_points = points;
initial_history = history;

% Original number of parameters
n_initial_parameters = size(history, 2);

% Optional initial conditions
if nargin < 3
    initial_conditions = zeros(n_initial_parameters, 1) - 1;
    step_size = 0.2;
elseif nargin == 3
    initial_conditions = varargin{1};
    step_size = 0.2;
elseif nargin == 4
    initial_conditions = varargin{1};
    step_size = varargin{2};
else
    error('Must be at least 2 inputs. Option inputs are initial conditions and step size')
end

% Input points and initial conditions should be a column vector. 
% If it is a row, transpose
if isrow(points)
    points = points';
end
if isrow(initial_conditions)
    initial_conditions = initial_conditions';
end
% Validate that they are now columns
if ~iscolumn(points)
    error('points input must be a column vector')
end
if ~iscolumn(initial_conditions)
    error('initial conditions must be a column vector')
end

% First remove all parameters that are trivially -Inf
% i.e., if this event has occured, then the current process will
% not have an event
inhibitory_parameters = find(points'*history == 0);
[inhibitory_bins, ~] = find(history(:, inhibitory_parameters) == 1);

history(:, inhibitory_parameters) = [];
history(inhibitory_bins, :) = [];
points(inhibitory_bins) = [];

% List of non-inhibitory indices for mapping the excitatory indices back to
% the original index
index_mapping = setdiff(1:n_initial_parameters, inhibitory_parameters);

% Now remove all parameters taht are trivially +Inf
% i.e., if this event has occured, the the current process will
% definitely have an event.

excitatory_parameters = find(points'*history == sum(history));
[excitatory_bins, ~] = find(history(:, excitatory_parameters) == 1);

history(:, excitatory_parameters) = [];
history(excitatory_bins, :) = [];
points(excitatory_bins) = [];

% Create a final list of inhibitory, excitatory and estimated indices to
% format final parameter
inhibitory_indices = inhibitory_parameters;
excitatory_indices = index_mapping(excitatory_parameters);
estimated_indices = setdiff(index_mapping, excitatory_indices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flag for determining recursion or model fit error or no problems
flag=0;
% convergence criteria
epsilon=0.00001;
min_step_size = 0.01;
% number of segments and parameters to estimate
n_bins = length(points);
n_parameters = size(history, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outputs
parameters = initial_conditions(estimated_indices);

% Gradiant (F) and Jacobian (J)
J = zeros(n_parameters);
F = zeros(n_parameters, 1);

% CIF
lambda = zeros(n_bins, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1 == 1
    for bin = 1:n_bins
        % Calculate CIF
        lambda(bin) = exp(sum(parameters' .* history(bin, :)));
        % If probability of a shift grows too large, model fails
        % Initial conditions should be less than 0 to accomodate
        if lambda(bin) > 1
            flag=1;
        end
    end
    if flag==1
        break
    end
    % Calculate gradiant (F) and Jacobian (J) for newtons method
    for parameter1 = 1:n_parameters
        F(parameter1) = -sum(points .* history(:, parameter1) + history(:, parameter1) .* (lambda .* ...
            (points - 1)) ./ (1 - lambda));
        for parameter2 = 1:n_parameters
            J(parameter1, parameter2) = sum(history(:, parameter1) .* history(:, parameter2) .* (points - 1) .* ...
                lambda ./ (1 - lambda) .^ 2);
        end
    end
    % Calculate change in parameters. Magnitude is reduced to 20% to
    % prevent parameters from over shooting and jumping outside the support 
    % of the model
    change_in_parameters = J \ F;
    parameters = parameters + change_in_parameters * step_size; 
    % If any of the parameters are below -20, set flag 2
    % Step size too large, or inhibitory parameter
    if min(parameters) < -10
        flag = 2;
        break;
    end
    % When convergence is reached, exit loop
    if max(abs(change_in_parameters)) < epsilon 
        break
    end
end

% flag=0: No errors, calculate output
if flag == 0
    for bin = 1:n_bins
        lambda(bin) = exp(sum(parameters' .* history(bin,:))); 
    end
    likelihood = sum((points .* log(lambda) + (1 - points) .* log(1 - lambda)));
    final_parameters = zeros(1, n_initial_parameters);
    final_parameters(estimated_indices) = parameters;
    final_parameters(inhibitory_indices) = -Inf;
    final_parameters(excitatory_indices) = Inf;
elseif step_size > min_step_size
    % flag set - try with smaller step size until below minimum
    [likelihood, final_parameters] = PSGC(initial_points, initial_history, initial_conditions, step_size / 2);
else
    error('Error: Model does not converge. Possible inhibitory or excitatory parameter')
end

