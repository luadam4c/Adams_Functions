function fieldsTable = all_fields (inStruct, varargin)
%% Get all field values and names of a structure that satisfies specific conditions in cell arrays
% Usage: fieldsTable = all_fields (inStruct, varargin)
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
%
% Outputs:
%       fieldsTable - a table with columns:
%                       Name        - the names of selected fields
%                       Value       - the values of selected fields
%                       IndexOrig   - the original index of selected fields
%                   specified as a column cell array
%
% Arguments:
%       inStruct    - a structure with fieldsTable
%                   must be a scalar structure
%       varargin    - 'ValueFunc': a function that takes the field value
%                           as the sole argument and a boolean as an output
%                   must be a function handle
%                   default == empty function handle
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/spike2Mat2Text.m

% File History:
% 2019-09-02 Created by Adam Lu
% TODO: Extract code for regexp match from all_files.m into a function
% TODO: Add regexp match
% 

%% Hard-coded parameters

%% Default values for optional arguments
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
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValueFunc', valueFuncDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));

% Read from the Input Parser
parse(iP, inStruct, varargin{:});
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
fieldValues = cellfun(@(x) inStruct.(x), fieldNames, 'UniformOutput', false);

% Return original indices
indOriginal = find(toSelect);

%% Output
% Place in a table
fieldsTable = table(fieldNames, fieldValues, indOriginal, ...
                    'VariableNames', {'Name', 'Value', 'IndexOrig'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%