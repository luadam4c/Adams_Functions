function indices = create_indices (endPoints, varargin)
%% Creates indices from endpoints (starting and ending indices)
% Usage: indices = create_indices (endPoints, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       indices     - indices for each pair of idxStart and idxEnd
%                   specified as a numeric vector 
%                       or a cell array of numeric vectors
% Arguments:
%       endPoints   - the starting and ending indices
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%       varargin    - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%
% Used by:
%       cd/extract_subvectors.m
%       cd/parse_pulse_response.m

% File History:
% 2018-12-17 Created by Adam Lu
% 2018-12-18 Added 'ForceCellOutput' as an optional argument
% 

%% Hard-coded parameters

%% Default values for optional arguments
forceCellOutputDefault = false;  % don't force output as a cell array by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'endPoints', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, endPoints, varargin{:});
forceCellOutput = iP.Results.ForceCellOutput;

%% Do the job
% Extract the starting and ending indices
[idxStart, idxEnd] = ...
    argfun(@(x) extract_elements(endPoints, x), 'first', 'last');

% Construct vectors of indices
if numel(idxStart) > 1 || numel(idxEnd) > 1
    indices = arrayfun(@(x, y) transpose(x:y), idxStart, idxEnd, ...
                    'UniformOutput', false);
else
    indices = transpose(idxStart:idxEnd);
end

% Force as cell array output if requested
if forceCellOutput && ~iscell(indices)
    indices = {indices};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%