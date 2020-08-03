function [yLimits, yTicks] = m3ha_decide_on_ylimits (measureTitle, varargin)
%% Decides on the y axis limits based on the measure title
% Usage: [yLimits, yTicks] = m3ha_decide_on_ylimits (measureTitle, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       yLimits     - y axis limits
%                   specified as a TODO
%
% Arguments:
%       measureTitle    - measure title
%                       must be a TODO
%       plotType    - 'PlotType': plot type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '2d1' - violin plots, Figures 2, 7 & 8
%                       '2d2' - violin plots, Figure 4
%                       '3d1' - bar plots, Figure 2
%                       '3d2' - bar plots, Figure 4
%                   default == '2d1'
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_plot_figure07.m
%       cd/m3ha_plot_grouped_scatter.m
%       cd/m3ha_plot_violin.m
%       cd/m3ha_plot_bar3.m

% File History:
% 2020-05-11 Moved from m3ha_plot_violin.m
% 2020-07-29 Added 'PlotType' as an optional argument

%% Hard-coded parameters
validPlotTypes = {'2d1', '2d2', '3d1', '3d2'};

%% Default values for optional arguments
plotTypeDefault  = '2d1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'plotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));

% Read from the Input Parser
parse(iP, varargin{:});
plotType = validatestring(iP.Results.plotType, validPlotTypes);

%% Preparation
% Initialize outputs
yLimits = [];
yTicks = [];

%% Do the job
switch measureTitle
    case {'LTS probability', 'Burst probability'}
        switch plotType
            case {'2d1', '2d2'}
                yLimits = [0, 1];
            case {'3d1', '3d2'}
                yLimits = [];
            otherwise
        end
    case {'LTS onset time (ms)', 'Burst onset time (ms)'}
        switch plotType
            case '2d1'
                yLimits = [0, 3500];
            case '2d2'
                yLimits = [0, 2000];
            case '3d1'
                yLimits = [0, 4000];
            case '3d2'
                yLimits = [0, 3000];
            otherwise
        end
    case {'LTS onset time (s)', 'Burst onset time (s)'}
        switch plotType
            case '2d1'
                yLimits = [0, 3.5];
            case '2d2'
                yLimits = [0, 2];
            case '3d1'
                yLimits = [0, 4];
            case '3d2'
                yLimits = [0, 3];
            otherwise
        end
    case 'Spikes Per LTS'
        switch plotType
            case '2d1'
                yLimits = [0, 12];
            case '2d2'
                yLimits = [0, 6.5];
            case '3d1'
                yLimits = [0, 8];
            case '3d2'
                yLimits = [0, 6];
            otherwise
        end
    case 'LTS maximum slope (V/s)'
        switch plotType
            case '2d1'
                yLimits = [0, 5];
            case '2d2'
                yLimits = [0, 8];
                yTicks = [0:2:8];
            case {'3d1', '3d2'}
                yLimits = [0, 4];
            otherwise
        end
    case 'LTS amplitude (mV)'
        switch plotType
            case '2d1'
                yLimits = [-75, -35];
            case '2d2'
                yLimits = [-75, -30];
                yTicks = -75:10:-35;
            case {'3d1', '3d2'}
                yLimits = [-75, -0];
            otherwise
        end
    case 'Oscillation Probability'
        yLimits = [0, 1];
    case 'Oscillation Period (sec)'
        % yLimits = [0.1, 1];
        yLimits = [0.3, 1.3];
    case 'Oscillatory Index'
        % yLimits = [0.1, 0.8];
        yLimits = [0.3, 1];
    case {'Active Cells (%)', 'Active TC Cells (%)'}
        yLimits = [0, 100];
    case 'Half Activation Time (sec)'
        % yLimits = [0, 10];
        yLimits = [0, 9];
    case 'Oscillation Duration (sec)'
        yLimits = [0, 30];
    otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%