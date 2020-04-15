function hold_off (wasHold)
%% Holds off based on previous status
% Usage: hold_off (wasHold)
% Explanation:
%       TODO
%
% Example(s):
%       wasHold = hold_on;
%       hold_off(wasHold);
%
% Arguments:
%       wasHold     - whether the current axes was held on
%                   must be numeric/logical 1 (true) or 0 (false)
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m
%       cd/plot_autocorrelogram.m
%       cd/plot_ball_stick.m
%       cd/plot_chevron.m
%       cd/plot_grouped_jitter.m
%       cd/plot_raw_multiunit.m
%       cd/plot_spike_density_multiunit.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-10-02 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'wasHold', ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, wasHold);

%% Do the job
if ~wasHold
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%