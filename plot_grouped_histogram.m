function plot_grouped_histogram(figname, stats, grouping, grouping_labels, xLabel, xUnits, titleStr, varargin)
%% Plot a grouped histogram
% Usage: plot_grouped_histogram(figname, stats, grouping, grouping_labels, xLabel, xUnits, titleStr, varargin)
%
% Requires:
%		/home/Matlab/Adams_Functions/histg.m
%
% Used by:
%		/media/adamX/Paula_IEIs/paula_iei4.m
%
% 20171211 - Created

%% Default values for optional arguments
yLabelDefault = 'Count';
xLimitsDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'plot_grouped_histogram';

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after 2016B
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric', 'categorical', ...
        'datetime', 'duration'}, {'vector', 'numel', 2}));

% Read from the Input Parser
parse(iP, varargin{:});
yLabel = iP.Results.YLabel;
xLimits = iP.Results.XLimits;

%% Plot and save histogram
h = figure('Visible', 'off');
clf(h);
histg(stats, grouping);
if ~isempty(xLimits)
    xlim(xLimits);
end
legend(grouping_labels, 'Interpreter', 'none', 'location', 'eastoutside');    
if ~isempty(xUnits)
    xlabel([xLabel, ' (', xUnits, ')']);
else
    xlabel(xLabel);
end
ylabel(yLabel);
title(titleStr, 'Interpreter', 'none');
saveas(h, figname, 'png');
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

