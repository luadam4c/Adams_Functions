function indZeros = find_zeros (vecs, varargin)
%% Find the indices where the vectors are closest to zero
% Usage: indZeros = find_zeros (vecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       indZeros    - TODO: Description of indZeros
%                   specified as a TODO
%
% Arguments:
%       vecs        -  TODO: Description of vecs
%                   must be a TODO
%       n           - (opt) Number of zeros to find
%                   must be a positive integer scalar
%                   default == [] (find all)
%       direction   - (opt) Search direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'first' - Search forward from the beginning
%                       'last'  - Search backward from the end
%                   default == 'first'
%       varargin    - 'ReturnNan': Return NaN instead of empty if nothing found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/set_default_flag.m
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m

% File History:
% 2020-04-15 Created by Adam Lu
% 2020-04-22 Added n and direction as optional arguments

%% Hard-coded parameters
validDirections = {'first', 'last'};

%% Default values for optional arguments
nDefault = [];              % default number of nonzero elements to find
directionDefault = 'first'; % default search direction
returnNanDefault = [];      % whether to return NaN instead of empty 
                            %   if nothing found by default

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
addRequired(iP, 'vecs');

% Add parameter-value pairs to the Input Parser
% Add optional inputs to the Input Parser
addOptional(iP, 'n', nDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['n must be either empty ', ...
                    'or a positive integer scalar!']));
addOptional(iP, 'direction', directionDefault, ...
    @(x) any(validatestring(x, validDirections)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ReturnNan', returnNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vecs, varargin{:});
n = iP.Results.n;
direction = validatestring(iP.Results.direction, validDirections);
returnNan = iP.Results.ReturnNan;

%% Preparation
% Force as column vectors
vecs = force_column_vector(vecs);

% Decide on uniformOutput
uniformOutput = set_default_flag([], ~isempty(n) && n == 1);

% Decide on returnNan
returnNan = set_default_flag(returnNan, uniformOutput);

%% Do the job
if iscell(vecs)
    indZeros = cellfun(@(x) find_zeros_helper(x, n, direction, returnNan), ...
                        vecs, 'UniformOutput', uniformOutput);
else
    indZeros = find_zeros_helper(vecs, n, direction, returnNan);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indZeros = find_zeros_helper(vec, n, direction, returnNan)

% Find the left and right of each consecutive pair of sample points
vecLeft = vec(1:end-1);
vecRight = vec(2:end);

% Find the indices with vec values closest to zero
indLeftZeros = find(vecLeft .* vecRight <= 0);

% Return if empty
if isempty(indLeftZeros)
    if returnNan
        indZeros = NaN;
    else
        indZeros = [];
    end
    return
end

% Get the next index
indRightZeros = indLeftZeros + 1;

% Choose the index closest to zero
indZeros = arrayfun(@(a, b) choose_closest_to_value(vec, a, b, 0), ...
                    indLeftZeros, indRightZeros);

% Restrict to the number of desired zeros
if ~isempty(n) && n < numel(indZeros)
    switch direction
        case 'first'
            indZeros = indZeros(1:n);
        case 'last'
            indZeros = indZeros((end-n+1):end);
        otherwise
            error('direction unrecognized!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = choose_closest_to_value (vec, idx1, idx2, value)

if abs(vec(idx1) - value) <= abs(vec(idx2) - value)
    idx = idx1;
else
    idx = idx2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%