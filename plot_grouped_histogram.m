function h = plot_grouped_histogram(figName, stats, grouping, grouping_labels, xLabel, xUnits, titleStr, varargin)
%% Plot a grouped histogram
% Usage: h = plot_grouped_histogram(figName, stats, grouping, grouping_labels, xLabel, xUnits, titleStr, varargin)
% Explanation:
%       TODO
%       Note: The bar() function is actually used for the main histogram
% Example(s):
%       TODO
% Outputs:
%       h           - the histogram returned as a Bar object
%                   specified as a Patch (R2015a) or Bar (R2017a) object
%
% Arguments: TODO
%       stats       - data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'OutFolder': directory to save figure, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == pwd
% Requires:
%       cd/create_error_for_nargin.m
%       cd/histg.m
%
% Used by:
% TODO:
%       cd/ZG_fit_IEI_distributions.m
%
% TODO:
%        /media/adamX/Paula_IEIs/paula_iei4.m
%       /home/Matlab/Marks_Functions/paula/Oct2017/freqsPostJustinPartTwo.m
%
% 2017-12-11 Created by Adam Lu
% 2018-05-18 Added outFolder as a parameter
% 2018-05-25 Now doesn't plot if stats is empty

%% Default values for optional arguments
yLabelDefault = 'Count';
xLimitsDefault = [];
outFolderDefault = '';          % default directory to save figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

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
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
yLabel = iP.Results.YLabel;
xLimits = iP.Results.XLimits;
outFolder = iP.Results.OutFolder;

% If the figure name is not a full path, create full path
if ~isempty(figName) && ~any(strfind(figName, filesep))
    % Set dependent argument defaults
    if isempty(outFolder)
        % Default output directory is present working directory
        outFolder = pwd;
    end

    % Create full figure file name
    figName = fullfile(outFolder, figName);
end

%% Plot and save histogram
h = figure('Visible', 'on');
clf(h);
if ~isempty(stats)
    histg(stats, grouping);
end
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
saveas(h, figName, 'png');
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

histg(stats, grouping);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
