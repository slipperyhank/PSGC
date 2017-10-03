

% In order to use the PSGC script, you require the following information
% about the EEG recording
%
% Inputs: 
%    data (matrix(float)) - EEG data to be analyzed. Can be either raw data
%       or ICA activations
%    boundary_markers (array(int)) - An array of indicies which mark the
%       beginning of a new segment / epoch. 
%    sampling_rate (int) - The sampling rate of the recording


% Here we consider a file that has already been preprocessed in eeglab. 
filename = 'P010_PPe_visual.set';
EEG = pop_loadset(filename);

% For this example, we are specifically interested in frontal and occipital
% independent components. At has been predetermined that components 6 and
% 17 are frontal, and components 1 and 5 are occipital. 
channels = [1, 5, 6, 17]; 
data = EEG.icaact(channels, :);
sampling_rate = EEG.srate;

% The PSGC program comes with a helper file the load the boundary_markers.
% The second input is the string name for boundary events. 
boundary_markers = eeglab_get_boundary_markers(EEG, 'boundary');


% The PSGC analysis also requires several parameters
%
% Parameters:
%   frequency band (array [lower, upper]): frequency band of interest
%   points_per_bin (int): number of data points to be grouped together in
%      one bin
%   bins_per_window (int): number of bins to be used in each history window
%   maximum_model_order (int): number of history windows to consider in
%      model selection



% Setting the frequency band.

% The PSGC analysis should be performed at a frequency band where there is
% an oscillator. Although it is difficult to know for sure, this can be
% tested by looking for small (or large) "bumps" in the spectral density of
% the EEG recordings. 

% A helper function is available to plot the estimated power spectral
% density for each channel in the data
explore_frequency_spectrum(data, sampling_rate);

% From the plots, there is a clear "bump" in the PSD. There is some
% variability in the bump from channel to channel, I would say they are [6,
% 10], [6, 13], [6, 10] and [5, 10].

% We want to select a frequency band that captures as much as possible, but
% also we want to keep it as narrow as possible. In this case, I will
% choose [6, 10]. Generally, I don't like having a frequency band more than
% 4-5 Hz in width. 
frequency_band = [6, 10];



% Setting the points_per_bin

% The points_per_bin is the temporal resolution for identifying phase phase
% shift events. The larger the value of points_per_bin, the faster the 
% analysis will run. However, for the PSGC model, it must be selected to be
% small enough such that the probability of observing two phase shifts in a
% bin is negligible. 

% Testing several values of points_per_bin
points_per_bin = 5;
point_process = make_point_process(data, boundary_markers, sampling_rate, frequency_band, points_per_bin);
max(max(point_process))

points_per_bin = 10;
point_process = make_point_process(data, boundary_markers, sampling_rate, frequency_band, points_per_bin);
max(max(point_process))

points_per_bin = 20;
point_process = make_point_process(data, boundary_markers, sampling_rate, frequency_band, points_per_bin);
max(max(point_process))

% At 20 points per bin, we see bins with two phase shift events, and thus
% 20 is too high. For 10 points per bin, we see a maximum of 1 shift per
% bin. We will set it to 10 for now, but if we are having model convergence
% issues, we may want to reduce it to 5 or lower. 
points_per_bin = 10;

% Setting bins_per_window and maximum model order. Bins per window and 
% maximum model order define the time scale that we expect to see one phase
% shift affecting another. 

% I typically choose bins_per_window so that a time window is roughly the
% length of one cycle of my frequency of interest. Here, we have a band
% centred at 8 Hz, and 1 cycle of 8 Hz takes 125 ms. With a sampling rate
% of 512, each point is approximately 2 ms and there are 10 points per bin,
% meaning each bin is approximate 20 ms. Setting bins_per_window = 6, we
% get that each window is approximately 120 ms, or one cycle of 8 Hz
% signel. 
bins_per_window = 6;


% To set the maximum model order, consider the total amount of time it
% should take to complete a task, or respond to a stimuli. If we expect the
% participant to respond / process the information within 500-600ms, we can
% select a maximum model order of 5, so that the history contains up to 5 *
% 120 = 600 ms. If instead we think it should take an entire second, we may
% consider m = 8, so that the history contains up to 8 * 120 = 960 ms. 
maximum_model_order = 8;



% Now that we have loaded the data and set the parameters, we can perform
% the PSGC analysis.
[pval_matrix, parameter_estimates, statistic_matrix, model_order] = psgc_analysis(data, boundary_markers, sampling_rate, frequency_band, points_per_bin, bins_per_window, maximum_model_order);

% There are four outputs from the PSGC analysis:
%
% pval_matrix (matrix(float)) (n_channels by n_channels): A matrix of
%    p-values, uncorrected for multiple comparisons. The p-value in the 
%    [i, j] entry of the matrix tests the hypothesis of no directed 
%    connection from channel j to channel i.
%
% parameter_estimates (cell array(matrix)): Parameter estimates for the
%    PSGC model for each individual channel. There are 1 + M * n_channel
%    parameters in the model, where M is the model order. The first
%    parameter is an intercept. The next M parameters represent the
%    influence the first channel has in each of the M time windows. The
%    second M values represent the influence the second channel has, etc. 
%
% statistic_matrix (matrix(float)): The chi-square values corresponding to
%    the pvalues in pval_matrix. Only used for post-hoc analysis. 
%
% model_order (array(n_channels)): The optimial model order for each
%    channel. Also the degrees of freedom for testing if another channel
%    has a directed influence over the channel being modelled.

% From the pval_matrix we see that pval_matrix[1, 3] = 0.0032, meaning that
% channel 3 has a significant effect on channel 1. Recall that channel 3 is
% a frontal channel and channel 1 is occipital. There are no other
% significant connections. Note, that if we perform a correction for the 8
% possible frontal to occipital connections we test, we get a pvalue of
% 0.0256 that is still significant. 

% We can further explore this significant interaction by looking at the
% model parameters for channel 3. % The influence of channel 1 on channel 3
% is described by parameters 2:8

parameter_estimates{3}(2:8)

% We see that channel 1 has an excitatory effect on channel 3, and that
% this effect is strongest in time windows 2 (0.18) and 3 (0.22). This
% means that when channel 1 experiences a phase shift, channel 3 is more
% likely to experience a phase shift 125 - 375 ms afterwards. 












