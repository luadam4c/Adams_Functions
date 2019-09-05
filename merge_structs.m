function structNew = merge_structs (struct1, struct2, varargin)
%% Merges two scalar structures
% Usage: structNew = merge_structs (struct1, struct2, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       structNew   - merged structure
%                   specified as a scalar structure
% Arguments:    
%       struct1     - first structure
%                   must be a scalar structure
%       struct2     - second structure
%                   must be a scalar structure
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:    
%       cd/fit_and_estimate_passive_params.m
%       cd/parse_stim.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting42.m

% File History:
% 2018-10-15 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'struct1', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addRequired(iP, 'struct2', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, struct1, struct2, varargin{:});

%% Do the job
% Initialize the new structure as the first structure
structNew = struct1;

% Add fields to this structure
allFields = fieldnames(struct2);
for iField = 1:numel(allFields)
    % Get this field name
    fieldNameThis = allFields{iField};

    % Copy over fields from the second structure
    structNew.(fieldNameThis) = struct2.(fieldNameThis);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%