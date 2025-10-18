function wasHold = hold_on(varargin)
%% Holds on and returns previous status
% Usage: wasHold = hold_on(axHandle (opt))
% Explanation:
%       TODO
%
% Example(s):
%       wasHold = hold_on;
%       hold_off(wasHold);
%
% Arguments:
%       axHandle    - (opt) axes handle to work with
%                   must be a empty or a axes object handle
%                   default == gca
%
% Outputs:
%       wasHold     - whether the current axes was held on
%                   specified as a logical scalar
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
%       cd/plot_test_result.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m

% File History:
% 2019-10-02 Created by Adam Lu
% 2025-09-17 Added axHandle as an optional argument

% Default values for optional arguments
axHandleDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add optional inputs to the Input Parser
addOptional(iP, 'axHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, varargin{:});
axHandle = iP.Results.axHandle;

%% Preparation
if isempty(axHandle)
    axHandle = gca;
end

%% Do the job
if ~ishold(axHandle)
    wasHold = false;
    hold(axHandle, 'on');
else
    wasHold = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%