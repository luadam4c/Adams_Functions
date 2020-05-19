function yLimits = m3ha_decide_on_ylimits (measureTitle)
%% Decides on the y axis limits based on the measure title
% Usage: yLimits = m3ha_decide_on_ylimits (measureTitle)
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
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_plot_figure07.m
%       cd/m3ha_plot_violin.m

% File History:
% 2020-05-11 Moved from m3ha_plot_violin.m
% 

%% Hard-coded parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
switch measureTitle
    case 'LTS probability'
        ylim([0, 1]);
    case 'LTS onset time (ms)'
        ylim([0, 2000]);
    case 'Spikes Per LTS'
        ylim([0, 6.5]);
    case 'LTS maximum slope (V/s)'
        % ylim([0, 5]);
        ylim([0, 8]);
        yticks([0:2:8]);
    case 'LTS amplitude (mV)'
        % ylim([-75, -45]);
        % yticks(-75:10:-45);
        ylim([-75, -30]);
        yticks(-75:10:-35);
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
        yLimits = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%