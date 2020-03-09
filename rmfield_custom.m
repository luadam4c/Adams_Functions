function structArray = rmfield_custom (structArray, fields, varargin)
%% Removes field(s) from a structure array only if the field(s) exists
% Usage: structArray = rmfield_custom (structArray, fields, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       rmfield_custom(myScalarStruct, ["a", "b"])
%       rmfield_custom(myScalarStruct, {'b', 'c'})
%       rmfield_custom(myScalarStruct, {'e', 'f'})
%
% Outputs:
%       structArray - returned structure array
%                   specified as a struct array
%
% Arguments:
%       structArray - structure array
%                   must be a struct array
%       fields      - fields to remove
%                   must be a character vector or a string array or 
%                       a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/all_subdirs.m
%       cd/array_fun.m
%       cd/compute_bins.m

% File History:
% 2019-11-13 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'structArray', ...
    @(x) validateattributes(x, {'struct'}, {'3d'}));
addRequired(iP, 'fields', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['fields must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, structArray, fields, varargin{:});
% param1 = iP.Results.param1;

%% Preparation
% Force character vectors as a cell array
if ischar(fields)
    fields = {fields};
end

%% Do the job
% Check whether each field is in the structure array
isField = isfield(structArray, fields);

% Remove only the fields in the structure array
fieldsToRemove = fields(isField);

% Remove each field
structArray = rmfield(structArray, fieldsToRemove);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%