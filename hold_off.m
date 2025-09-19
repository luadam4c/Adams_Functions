function hold_off (wasHold, varargin)
%% Holds off based on previous status
% Usage: hold_off (wasHold, axHandle (opt))
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
%       axHandle    - (opt) axes handle to work with
%                   must be a empty or a axes object handle
%                   default == gca
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compare_events_pre_post_stim.m
%       cd/create_plot_movie.m
%       cd/plot_autocorrelogram.m
%       cd/plot_ball_stick.m
%       cd/plot_chevron.m
%       cd/plot_correlation_coefficient.m
%       cd/plot_grouped_histogram.m
%       cd/plot_grouped_jitter.m
%       cd/plot_grouped_scatter.m
%       cd/plot_selected.m
%       cd/plot_spectrogram.m
%       cd/plot_spike_density_multiunit.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m

% File History:
% 2019-10-02 Created by Adam Lu
% 2025-09-17 Added axHandle as an optional argument

% Default values for optional arguments
axHandleDefault = gca;

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

% Add optional inputs to the Input Parser
addOptional(iP, 'axHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, wasHold, varargin{:});
axHandle = iP.Results.axHandle;

%% Preparation
if ~isempty(axHandle)
    axes(axHandle);
end

%% Do the job
if ~wasHold
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%