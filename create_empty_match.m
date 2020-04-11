function empty = create_empty_match (array, varargin)
%% Creates an empty array that matches a given array
% Usage: empty = create_empty_match (array, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       empty       - empty array matched to input
%
% Arguments:
%       array       - input array to match
%       varargin    - 'NRows': number of rows
%                   must be a nonnegative integer scalar
%                   default == size(array, 1)
%                   - 'NColumns': number of columns
%                   must be a nonnegative integer scalar
%                   default == size(array, 2)
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/all_ordered_pairs.m
%       cd/compute_combined_trace.m
%       cd/extract_subvectors.m
%       cd/force_matrix.m

% File History:
% 2019-01-03 Moved from extract_subvectors.m
% 2019-08-21 Defaults duration arrays to NaN minutes
% 2019-08-21 Defaults graphics handle arrays to gobjects
% 2019-10-03 Now creates a cell array of empty strings to match cellstrs
% 2020-02-04 Fixed the case for structures
% 

%% Hard-coded parameters

%% Default values for optional arguments
nRowsDefault = [];          % set later
nColumnsDefault = [];       % set later

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
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NRows', nRowsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));
addParameter(iP, 'NColumns', nColumnsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));

% Read from the Input Parser
parse(iP, array, varargin{:});
nRows = iP.Results.NRows;
nColumns = iP.Results.NColumns;

%% Preparation
% If not provided, get dimensions
if isempty(nRows)
    nRows = size(array, 1);
end
if isempty(nColumns)
    nColumns = size(array, 2);
end

%% Do the job
% Construct the empty array according to type
if isnumeric(array)
    empty = NaN(nRows, nColumns);
elseif isduration(array)
    empty = minutes(NaN(nRows, nColumns));
elseif islogical(array)
    empty = false(nRows, nColumns);
elseif iscellstr(array)
    empty = arrayfun(@(x) '', ones(nRows, nColumns), 'UniformOutput', false);
elseif iscell(array)
    empty = cell(nRows, nColumns);
elseif isstruct(array)
    empty = repmat(array([]), [nRows, nColumns]);
elseif isdatetime(array)
    empty = NaT(nRows, nColumns);
elseif isgraphics(array)
    empty = gobjects(nRows, nColumns);
elseif isa(array, 'table')
    % TODO: Not tested
    empty = table.empty(nRows, nColumns);
elseif isa(array, 'function_handle')
    % TODO: Not tested
    empty = function_handle.empty(nRows, nColumns);
else
    error('Not implemented yet!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%