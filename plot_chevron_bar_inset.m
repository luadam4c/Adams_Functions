function handles = plot_chevron_bar_inset (data, varargin)
%% Plots a bar plot comparing against the first group
% Usage: handles = plot_chevron_bar_inset (data, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles - handles structure with fields:
%                   bars    - bar objects (bars for each group or 
%                                           column is one bar object)
%                           specified as a vector of Bar object handles
%                   lines   - error bar lines
%                               1st dim: connecting (1), upper limit (2), lower limit (3)
%                               2nd dim: sample number
%                               3rd dim: group number
%                           specified as an array of Primitive Line object handles
%                   fig     - figure handle
%                           specified as a Figure object handle
%               specified as a structure
%
% Arguments:
%       data        - data table or data vectors
%                   Note: The dimension with fewer elements is taken as 
%                           the parameter
%                   must be a table or a numeric array
%                       or a cell array of numeric vectors
%       varargin    - 'PLimits': limits of parameter axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [0.5, nGroups + 0.5]
%                   - 'ReadoutLimits': limits of readout axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [0, Inf]
%                   - 'PTickLabels': tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == group labels detected from table
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == 'suppress'
%                   - 'ReadoutLabel': label for the readout
%                   must be a string scalar or a character vector
%                   default == 'suppress'
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == '% Baseline'
%                   - Any other parameter-value pair for plot_bar()
%
% Requires:
%       cd/compute_stats.m
%       cd/create_error_for_nargin.m
%       cd/force_data_as_matrix.m
%       cd/force_row_vector.m
%       cd/plot_bar.m
%
% Used by:
%       cd/m3ha_oscillations_analyze.m
%       cd/plot_measures.m

% File History:
% 2020-02-21 Moved from m3ha_oscillations_analyze.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
pLimitsDefault = [];                % set later
readoutLimitsDefault = [0, Inf];
pTickLabelsDefault = {};            % set later
pLabelDefault = 'suppress';
readoutLabelDefault = 'suppress';
figTitleDefault = '% Baseline';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'data', ...
    @(x) validateattributes(x, {'numeric', 'cell', 'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PLimits', pLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'ReadoutLimits', readoutLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, data, varargin{:});
pLimits = iP.Results.PLimits;
readoutLimits = iP.Results.ReadoutLimits;
pTickLabels = iP.Results.PTickLabels;
pLabel = iP.Results.PLabel;
readoutLabel = iP.Results.ReadoutLabel;
figTitle = iP.Results.FigTitle;

% Keep unmatched arguments for the plot_bar() function
otherArguments = iP.Unmatched;

%% Preparation
% Force data values as a numeric matrix 
%   where each group is a column and each row is a sample
[dataValues, groupLabelsAuto, sampleLabelsAuto] = force_data_as_matrix(data);

% Set default p tick labels
if isempty(pTickLabels)
    pTickLabels = groupLabelsAuto;
end

% Set default p limits
if isempty(pLimits)
    nGroups = size(dataValues, 2);
    pLimits = [0.5, nGroups + 0.5];
end

%% Compute
% Compute means
means = compute_stats(dataValues, 'mean', 'IgnoreNan', true);

% Compute means normalized against the first group
normalizedMeans = (means ./ means(1)) .* 100;

% Force as a row so that bars are different groups
normalizedMeans = force_row_vector(normalizedMeans);

%% Plot
% Plot bar inset
handles = plot_bar(normalizedMeans, 'PLimits', pLimits, ...
                    'ReadoutLimits', readoutLimits, ...
                    'PTickLabels', pTickLabels, ...
                    'ReadoutLabel', readoutLabel, 'PLabel', pLabel, ...
                    'FigTitle', figTitle, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%