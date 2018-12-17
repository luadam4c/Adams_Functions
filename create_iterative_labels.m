function labels = create_iterative_labels (nLabels, varargin)
%% Creates iterative labels (1 through nLabels) with optional prefix or suffix
% Usage: labels = create_iterative_labels (nLabels, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       labels     - labels created
%                   specified as a cell array of character vectors
% Arguments:
%       nLabels     - number of labels
%                   must be a positive integer scalar
%       varargin    - 'Prefix': string to place before each number
%                   must be a TODO
%                   default == TODO
%                   - 'Suffix': string to place after each number
%                   must be a TODO
%                   default == TODO
%
% Used by:
%       cd/argfun.m

% File History:
% 2018-12-17 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
prefixDefault = [];     % no string to place before each number by default
suffixDefault = [];     % no string to place after each number by default

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
addRequired(iP, 'nLabels', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, nLabels, varargin{:});
prefix = iP.Results.Prefix;
suffix = iP.Results.Suffix;

%% Do the job
labels = arrayfun(@(x) [prefix, num2str(x), suffix], 1:nLabels, ...
                        'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%