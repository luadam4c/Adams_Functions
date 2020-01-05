function [fieldValue, fieldName] = ...
                first_matching_field (genStructs, candNames, varargin)
%% Extracts the first matching field/variable/property of a structure/table/property from a list of candidate names
% Usage: [fieldValue, fieldName] = ...
%               first_matching_field (genStructs, candNames, varargin)
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
%                       if not found, [] is returned
%       fieldName   - original field name
%                       if not found, '' is returned
%                   specified as a character array
%
% Arguments:
%       genStructs  - structures in the general sense
%                       (structures/tables/objects) to extract from
%                   must be a struct/table/object array
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
%       cd/is_field.m
%
% Used by:
%       cd/decide_on_geom_params.m
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
addRequired(iP, 'genStructs');
addRequired(iP, 'candNames', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['candNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, genStructs, candNames, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if ischar(candNames)
    fieldName = candNames;
    fieldValue = genStructs.(candNames);
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
        if is_field(genStructs, candNameThis)
            fieldName = candNameThis;
            fieldValue = genStructs.(candNameThis);
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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%