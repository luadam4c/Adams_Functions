function [dataAvg, groups] = m3ha_average_by_group (dataOrig, grouping, varargin)
%% Average data according to a grouping vector
% Usage: [dataAvg, groups] = m3ha_average_by_group (dataOrig, grouping, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       dataAvg     - data averaged
%                   specified as a numeric array
%       groups      - unique grouping values
%                   specified as a numeric vector
% Arguments:
%       dataOrig    - original data
%                   must be a numeric array or a cell array of numeric arrays
%       grouping   - holding voltages conditions
%                   must be a numeric vector
%       varargin    - 'VecNumberToAverage': vector number to average over trace
%                   must be a positive integer scalar
%                   default == 2
%
% Requires:
%       cd/compute_average_trace.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/extract_columns.m
%       cd/isnum.m
%
% Used by:
%       cd/m3ha_import_raw_traces.m

% File History:
% 2019-01-12 Adapted from code in m3ha_import_raw_traces.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
vecNumberToAverageDefault = 2;

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
addRequired(iP, 'dataOrig', ...
    @(x) assert(isnum(x) || iscell(x), ...
                ['dataOrig must be either a numeric array ', ...
                    'or a cell array!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VecNumberToAverage', vecNumberToAverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, dataOrig, varargin{:});
vecNumberToAverage = iP.Results.VecNumberToAverage;

%% Preparation
% Count the number of vectors
nVectors = count_vectors(dataOrig);

%% Do the job
% Separate each vector from the matrix
%   Note: dataSeparated will be a cell array of cell arrays of column vectors
dataSeparated = extract_columns(dataOrig, 1:nVectors, 'OutputMode', 'single');

% Find unique grouping values
groups = unique(grouping, 'sorted');

% Count the number of groups
nGroups = length(groups);

% Extract the vector of interest

dataSeparated{vecNumberToAverage}


% Group the voltage traces by unique grouping values, 
%   then average the grouped traces
% TODO: Use compute_average_trace with modifications to pass group
vVecsAveraged = cell(nGroups, 1);
for iVhold = 1:nGroups
    % Get the current grouping value
    vnow = groups(iVhold);

    % Collect all cpr traces with this grouping value
    vVecsGroupedThisVhold = vVecs(grouping == vnow);

    % Average the voltage traces from this grouping group
    vVecAveragedThis = compute_average_trace(vVecsGroupedThisVhold);

    % Save in arrays
    vVecsAveraged{iVhold} = vVecAveragedThis;
end
vVecs = vVecsAveraged;

% For time, current and conductance, 
%   take the first copy and repeat nGroups times
%   Note: Since pulse widths are not necessarily the same, 
%           current vectors should not be averaged
[tVecs, iVecs, gVecs] = ...
    argfun(@(x) repmat(x(1), nGroups, 1), tVecs, iVecs, gVecs);

% Re-combine the vectors for output
dataAvg = cellfun(@(x, y, z, w) horzcat(x, y, z, w), ...
                    tVecs, vVecs, iVecs, gVecs, ...
                    'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[tVecs, vVecs, iVecs, gVecs] = extract_columns(dataOrig, 1:4);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%