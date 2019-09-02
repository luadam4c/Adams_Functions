function [fieldValues, fieldNames] = all_fields (myStruct, varargin)
%% Get all field values and names of a structure that satisfies specific conditions in cell arrays
% Usage: [fieldValues, fieldNames] = all_fields (myStruct, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       [v, n] = all_fields(blab)
%       [v, n] = all_fields(myStruct)
%       [v, n] = all_fields(myStruct, 'ValueFunc', @isnumeric)
%       [v, n] = all_fields(myStruct, 'ValueFunc', @islogical)
%       [v, n] = all_fields(myStruct, 'ValueFunc', @isnum)
%       [v, n] = all_fields(myStruct, 'ValueFunc', @istext)
%       [v, n] = all_fields(myScalarStruct)
%
% Outputs:
%       fieldValues - the field values selected
%                   specified as a column cell array
%       fieldNames  - the field names selected
%                   specified as a column cell array
%
% Arguments:
%       myStruct    - a structure with fieldValues
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
addRequired(iP, 'myStruct', ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValueFunc', valueFuncDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));

% Read from the Input Parser
parse(iP, myStruct, varargin{:});
valueFunc = iP.Results.ValueFunc;

%% Do the job
% Get all the field names
allFieldNames = fieldnames(myStruct);

% Decide on whether to select the field based on the value
if ~isempty(valueFunc)
    toSelect = cellfun(@(x) valueFunc(myStruct.(x)), allFieldNames);
else
    toSelect = true(size(allFieldNames));
end

% Restrict the field names
fieldNames = allFieldNames(toSelect);

% Get all the field values
fieldValues = cellfun(@(x) myStruct.(x), fieldNames, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%