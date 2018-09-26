function h = plot_fields (structArray, varargin)
%% Plot all fields from a structure array as tuning curves
% Usage: h = plot_fields (structArray, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - figure handle for the created figure
%                   specified as a figure handle
% Arguments:    
%       structArray - a structure array containing scalar fields
%                   must be a 2-D structure array
%                   - 'PisLog': whether parameter values are to be plotted 
%                               log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == [false, false];
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'YLimits': limits of y axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == []
%                   - 'PTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == 'Parameter'
%                   - 'FieldLabels': label for the field
%                   must be a cell array of character vectors/strings
%                   default == field name
%                   - 'SingleColor': color when colsToPlot == 1
%                   must be a 3-element vector
%                   - 'FigTitles': titles for each figure
%                   must be a cell array of character vectors/strings
%                   default == [fieldLabel, ' vs. ', xLabel]
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       cd/plot_tuning_curve.m
%       cd/isfigtype.m
%
% Used by:    
%       cd/plot_all_abfs_dir.m

% File History:
% 2018-09-26 Created by Adam Lu
% TODO: Add 'XTickLabel' and 'XLabel' as optional arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
xislogDefault = [false, false];
xlimitsDefault = [];
ylimitsDefault = [];
xTicksDefault = [];
xTickLabelsDefault = {};
xLabelDefault = 'Parameter';
fieldLabelsDefault = {};
singleColorDefault = [0, 0, 1];
figTitlesDefault = {};          % set later
figNumberDefault = [];          % invisible figure by default
figNameDefault = '';
figTypesDefault = 'png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'structArray', ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XisLog', xislogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'XLimits', xlimitsDefault, ...
    @(x) ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', ylimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                            {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'XTicks', xTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'XTickLabels', xTickLabelsDefault, ...
    @(x) iscellstr(x) || isstring(x));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FieldLabels', fieldLabelsDefault, ...
    @(x) iscellstr(x) || isstring(x));
addParameter(iP, 'SingleColor', singleColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
addParameter(iP, 'FigTitles', figTitlesDefault, ...
    @(x) iscellstr(x) || isstring(x));
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, structArray, varargin{:});
xIsLog = iP.Results.XisLog;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xTicks = iP.Results.XTicks;
xTickLabels = iP.Results.XTickLabels;
xLabel = iP.Results.XLabel;
fieldLabels = iP.Results.FieldLabels;
singlecolor = iP.Results.SingleColor;
figTitles = iP.Results.FigTitles;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Check relationships between arguments
if numel(xTicks) ~= numel(xTickLabels)
    fprintf('PTicks and PTickLabels must have the same number of elements!\n');
    h = [];
    return
end

%% Preparation
% Count the number of entries
nEntries = length(structArray);

% Return if there are no entries
if nEntries == 0
    h = [];
    return;
end

% Create an vector for the parameter values
xValues = 1:nEntries;
if isempty(xTicks)
    xTicks = xValues;
end
if isempty(xTickLabels)
    xTickLabels = arrayfun(@(x) num2str(x), xValues, 'UniformOutput', false);
end

% Get all the fields of the structArray as a cell array
allFields = fieldnames(structArray);

% Count the number of fields
nFieldsOrig = numel(allFields);

% Take only the fields of the structArray that are numeric scalars
scalarStructArray = structArray;
for iField = 1:nFieldsOrig
    % Get the field name
    thisFieldName = allFields{iField};

    % Get the first instance of this field value
    thisFieldValue = structArray(1).(thisFieldName);

    % Remove this field if not a numeric scalar
    if ~isscalar(thisFieldValue) || ~isnumeric(thisFieldValue)
        fprintf(['Warning: the field %s is not a ', ...
                    'numeric scalar so will be removed!!\n'], ...
                    thisFieldName);
        scalarStructArray = rmfield(scalarStructArray, thisFieldName);
    end
end

% Get all the fields of the scalarStructArray as a cell array
allScalarFields = fieldnames(scalarStructArray);

% Count the number of fields
nFields = numel(allScalarFields);

% Return if there are no more fields
if nFields == 0
    h = [];
    return;
end

% Convert the data to a homogeneous array, with each column being a field
fieldData = table2array(struct2table(scalarStructArray));

%% Plot all fields
for iField = 1:nFields
    % Get the field value vector for this field
    field = fieldData(:, iField);

    % Set the field label for this field
    if ~isempty(fieldLabels)
        % Use the user-provided field label
        fieldLabel = fieldLabels{iField};
    else
        % Use the field name
        fieldLabel = allScalarFields{iField};
    end
        
    % Set the figure title
    if ~isempty(figTitles)
        % Use the user-provided figure title
        figTitle = figTitles{iField};
    else
        % Use the default
        if ~strcmpi(fieldLabel, 'suppress') && ~strcmpi(xLabel, 'suppress')
            figTitle = strrep([fieldLabel, ' vs. ', xLabel], '_', '\_');
        elseif ~strcmpi(fieldLabel, 'suppress')
            figTitle = strrep([fieldLabel, ' vs. parameter'], '_', '\_');
        else
            figTitle = 'Readout vs. parameter';
        end
    end
    
    % Create a figure
    h = figure;
    
    % Plot the tuning curve
    h = plot_tuning_curve(xValues, field, 'PisLog', xIsLog, ...
                        'XLimits', xLimits, 'YLimits', yLimits, ...
                        'PTicks', xValues, 'PTickLabels', xTickLabels, ...
                        'PLabel', 'suppress', ...
                        'ReadoutLabel', fieldLabel, ...
                        'SingleColor', singlecolor, ...
                        'FigTitle', figTitle, 'FigNumber', figNumber, ...
                        'FigName', figName, 'FigTypes', figtypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

