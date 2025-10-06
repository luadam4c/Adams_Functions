function [newStruct, colorMapStruct] = convert_colors_to_rgb (oldStruct, varargin)
%% Converts color name strings in a structure to RGB vectors.
% Usage: [newStruct, colorMapStruct] = convert_colors_to_rgb(oldStruct, varargin)
%
% Explanation:
%       This function iterates through the fields of an input structure. 
%       It identifies fields whose names contain a specific substring (e.g., 'color')
%       and converts their string values (e.g., 'red', 'darkblue') into 
%       [R, G, B] numeric vectors based on a provided or default color map.
%
% Example(s):
%       params.Plotting.colorPb = 'red';
%       params.Plotting.colorRr = 'blue';
%       params.Plotting.lineWidth = 2;
%       params.Analysis.someOtherParam = 'test';
%       newParams = convert_colors_to_rgb(params);
%       disp(newParams.Plotting);
%
% Outputs:
%       newStruct       - The modified structure with RGB vectors for color fields.
%                       specified as a structure
%       colorMapStruct  - The color map that was used for the conversion.
%                       specified as a structure
%
% Arguments:
%       oldStruct   - The input structure to process.
%                   must be a structure
%       varargin    - 'ColorSubStr': Substring to identify color fields.
%                   must be a string scalar or a character vector
%                   default == 'color'
%                   - 'ColorMapStruct': A structure mapping color names to RGB vectors.
%                   must be a structure
%                   default == (see default values below)
%
% Requires:
%       cd/decide_on_colormap.m 
%
% Used by:
%       \Shared\Code\vIRt-Moore\virt_moore.m
%       \Shared\Code\vIRt-Moore\virt_moore_monte_carlo.m

% File History:
% 2025-10-03 Created by Gemini from virt_moore.m

%% Default values for optional arguments
colorSubStrDefault = 'color';       % Default substring for identifying color fields
colorMapStructDefault = struct( ... % Default custom color map
    'red', [1, 0, 0], 'green', [0, 1, 0], 'blue', [0, 0, 1], ...
    'cyan', [0, 1, 1], 'magenta', [1, 0, 1], 'yellow', [1, 1, 0], ...
    'black', [0, 0, 0], 'white', [1, 1, 1], 'gray', [0.5, 0.5, 0.5], ...
    'lightgray', [0.85 0.85 0.85], ...
    'darkgreen', [0, 0.5, 0], 'darkblue', [0, 0, 0.5], ...
    'orange', [1, 0.65, 0], 'lightorange', [1.0, 0.8, 0.5], ...
    'pink', [1, 0.75, 0.8], 'aquamarine', [0.5, 1.0, 0.8], ...
    'default', [] ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;

% Add required inputs to the Input Parser
addRequired(iP, 'oldStruct', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ColorSubStr', colorSubStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColorMapStruct', colorMapStructDefault, ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Read from the Input Parser
parse(iP, oldStruct, varargin{:});
colorSubStr = iP.Results.ColorSubStr;
colorMapStruct = iP.Results.ColorMapStruct;

%% Do the job
% Create a copy of the input structure to modify
newStruct = oldStruct;

% Get all field names from the structure
allFields = fieldnames(newStruct);

% Note: A for-loop is used here because the decision to convert a field's 
%       value depends on its name. The `structfun` function applies a 
%       function to each field's value but does not provide the field's 
%       name as an input, making it unsuitable for this specific task.

% Iterate over each field
for iField = 1:numel(allFields)
    fieldName = allFields{iField};
    
    % Check if the field name contains the color substring (case-insensitive)
    if contains(fieldName, colorSubStr, 'IgnoreCase', true)
        originalValue = newStruct.(fieldName);
        
        % Check if the value is a string and not already an RGB vector
        if ischar(originalValue) || isstring(originalValue)
            % Convert the color string to an RGB vector
            newStruct.(fieldName) = get_color_rgb(originalValue, colorMapStruct);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rgb = get_color_rgb(colorStr, colorMap)
%% Converts a color string to an RGB vector using a map or a helper function.

% Convert the input string to a character vector for consistency
colorStr = char(colorStr);

% First, check if the color exists in the custom color map structure
if isfield(colorMap, colorStr)
    % Use custom color map
    rgb = colorMap.(colorStr);
elseif isempty(colorStr)
    % If the string is empty, return an empty value
    rgb = [];
else
    % Otherwise, use decide_on_colormap.m to determine the color
    rgb = decide_on_colormap(colorStr, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%