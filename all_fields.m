function output = all_fields (inStruct, varargin)
%% Get all field values and names of a structure that satisfies specific conditions in cell arrays
% Usage: output = all_fields (inStruct, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       all_fields(blab)
%       all_fields(myStruct)
%       all_fields(myStruct, 'ValueFunc', @isnumeric)
%       all_fields(myStruct, 'ValueFunc', @islogical)
%       all_fields(myStruct, 'ValueFunc', @isnum)
%       all_fields(myStruct, 'ValueFunc', @istext)
%       all_fields(myScalarStruct)
%       all_fields(myScalarStruct, 'OutputType', 'name')
%       all_fields(myScalarStruct, 'OutputType', 'value')
%       all_fields(myScalarStruct, 'OutputType', 'index')
%
% Outputs:
%       output      - a table with columns:
%                       Name        - the names of selected fields
%                       Value       - the values of selected fields
%                       IndexOrig   - the original index of selected fields
%                       or just the column requested
%                   specified as a column cell array
%
% Arguments:
%       inStruct    - a structure with fields
%                   must be a scalar structure
%       varargin    - 'OutputType': type of output
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'table'     - everything in a table
%                       'name'      - just the names in a cell array
%                       'value'     - just the values in a cell array
%                       'index'     - just the original indices
%                   default == 'table'
%                   - 'ValueFunc': a function that takes the field value
%                           as the sole argument and a boolean as an output
%                   must be a function handle
%                   default == empty function handle
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/parse_spike2_mat.m

% File History:
% 2019-09-02 Created by Adam Lu
% 2019-09-22 Added 'OutputType' as an optional argument
% 2019-09-29 Now accepts matfiles
% TODO: Extract code for regexp match from all_files.m into a function
% TODO: Add regexp match
% 

%% Hard-coded parameters
validOutputTypes = {'table', 'name', 'value', 'index'};

%% Default values for optional arguments
outputTypeDefault = 'table';
valueFuncDefault = function_handle.empty;

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
addRequired(iP, 'inStruct', ...
    @(x) validateattributes(x, {'struct', 'matlab.io.MatFile'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputType', outputTypeDefault, ...
    @(x) any(validatestring(x, validOutputTypes)));
addParameter(iP, 'ValueFunc', valueFuncDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));

% Read from the Input Parser
parse(iP, inStruct, varargin{:});
outputType = validatestring(iP.Results.OutputType, validOutputTypes);
valueFunc = iP.Results.ValueFunc;

%% Do the job
% Get all the field names
allFieldNames = fieldnames(inStruct);

% Decide on whether to select the field based on the value
if ~isempty(valueFunc)
    toSelect = cellfun(@(x) valueFunc(inStruct.(x)), allFieldNames);
else
    toSelect = true(size(allFieldNames));
end

% Restrict the field names
fieldNames = allFieldNames(toSelect);

% Get all the field values
if ~strcmp(outputType, 'name')
    fieldValues = cellfun(@(x) inStruct.(x), ...
                            fieldNames, 'UniformOutput', false);
end

% Return original indices
if ~strcmp(outputType, 'name') && ~strcmp(outputType, 'value')
    indOriginal = find(toSelect);
end

%% Output
switch outputType
    case 'table'
        % Place in a table
        output = table(fieldNames, fieldValues, indOriginal, ...
                            'VariableNames', {'Name', 'Value', 'IndexOrig'});
    case 'name'
        output = fieldNames;
    case 'value'
        output = fieldValues;
    case 'index'
        output = indOriginal;
    otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%