function [fieldValue, fieldName] = ...
                first_matching_field (structsOrtable, candNames, varargin)
%% Extracts the first matching field (or variable) of a structure (or table) from a list of candidate field (variable) names
% Usage: [fieldValue, fieldName] = ...
%               first_matching_field (structsOrtable, candNames, varargin)
% Explanation:
%       TODO
%
%       See also:
%           cd/extract_fields.m
%
% Example(s):
%       load_examples;
%       [fieldValue, fieldName] = first_matching_field(blab, {'Students', 'students'})
%       [fieldValue, fieldName] = first_matching_field(blab, {'None', 'Nil'})
%
% Outputs:
%       fieldValue  - extracted field value
%       fieldName   - original field name
%                   specified as a character array
%
% Arguments:
%       structsOrtable     - structures to extract from
%                   must be a struct array or a table
%       candNames   - candidate field names
%                   must be a character vector, a string array 
%                       or a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/is_var_in_table.m
%
% Used by:
%       cd/m3ha_compute_statistics.m
%       cd/plot_ball_stick.m

% File History:
% 2019-12-26 Created by Adam Lu
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
addRequired(iP, 'structsOrtable', ...
    @(x) validateattributes(x, {'struct', 'table'}, {'2d'}));
addRequired(iP, 'candNames', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['candNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, structsOrtable, candNames, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if ischar(candNames)
    fieldName = candNames;
    fieldValue = structsOrtable.(candNames);
else
    % Count the number of candidates
    nCands = numel(candNames);

    % Test each candidate field name in turn
    found = false;
    for iCand = 1:nCands
        % Extract the name
        if iscell(candNames)
            candNameThis = candNames{iCand};
        else
            candNameThis = candNames(iCand);
        end

        % Test whether it's a valid field
        if is_field_or_var(structsOrtable, candNameThis)
            fieldName = candNameThis;
            fieldValue = structsOrtable.(candNameThis);
            found = true;
            break;
        end
    end

    % If nothing found, return empty values
    if ~found
        fieldName = '';
        fieldValue = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isFieldOrVar = is_field_or_var(structsOrtable, candName)

if isstruct(structsOrtable)
    isFieldOrVar = isfield(structsOrtable, candName);
else
    isFieldOrVar = is_var_in_table(candName, structsOrtable);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%