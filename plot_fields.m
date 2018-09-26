function h = plot_fields (structArray, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
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
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/plot_tuning_curve.m
%
% Used by:    
%       /TODO:dir/TODO:file

% File History:
% 2018-09-26 Created by Adam Lu
% TODO: Add 'XTickLabel' and 'XLabel' as optional arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
% TODO

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

% Read from the Input Parser
parse(iP, structArray, varargin{:});
% TODO

% Check relationships between arguments
% TODO

%% Preparation
% TODO

% Count the number of entries
nEntries = length(structArray);

% Return if there are no entries
if nEntries == 0
    h = [];
    return;
end

% Create an vector for the parameter values
pValues = 1:nEntries;
pTickLabels = arrayfun(@(x) num2str(x), pValues, 'UniformOutput', false);

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
    % Get the readout vector for this field
    readout = fieldData(:, iField);

    % Get the readout label for this field
    readoutLabel = allScalarFields{iField};
    
    % Create a figure
    h = figure;
    
    % Plot the tuning curve
    h = plot_tuning_curve(pValues, readout, ...
                        'PTicks', pValues, ...
                        'PTickLabels', pTickLabels, ...
                        'PLabel', 'suppress', ...
                        'ReadoutLabel', readoutLabel);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

