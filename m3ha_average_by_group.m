function [dataAvg, groups] = m3ha_average_by_group (dataOrig, grouping, varargin)
%% Average data according to a grouping vector
% Usage: [dataAvg, groups] = m3ha_average_by_group (dataOrig, grouping, varargin)
% Explanation:
%       TODO
% Example(s):
%       m3ha_average_by_group({magic(3), magic(3) + 1, magic(3) + 2}, {'a', 'b', 'b'}, 'ColNumToAverage', 2:3)
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
%       varargin    - 'ColNumToAverage': column number(s) to average
%                   must be empty or a positive integer vector
%                   default == transpose(1:min(nColumns))
%
% Requires:
%       cd/compute_average_trace.m
%       cd/compute_combined_trace.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/extract_columns.m
%       cd/isnum.m
%
% Used by:
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-01-12 Adapted from code in m3ha_import_raw_traces.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
colNumToAverageDefault = [];

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
addParameter(iP, 'ColNumToAverage', colNumToAverageDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x), ...
                ['ColNumToAverage must be either empty or , ', ...
                    'a positive integer vector!']));

% Read from the Input Parser
parse(iP, dataOrig, varargin{:});
colNumToAverage = iP.Results.ColNumToAverage;

%% Preparation
% Count the number of columns for each sweep
nColumns = count_vectors(dataOrig);

% Generate a vector of column numbers
allColNums = create_indices('IndexEnd', min(nColumns));

% Average all columns by default
if isempty(colNumToAverage)
    colNumToAverage = allColNums;
end

% Compute other column numbers
colNumOther = setdiff(allColNums, colNumToAverage);

% Find unique grouping values
groups = unique(grouping, 'sorted');

% Count the number of groups
nGroups = length(groups);

% Generate a vector of group numbers
allGroupNums = create_indices('IndexEnd', nGroups);

%% Do the job
% Extract each column from all sweeps separately
%   Note: dataSeparated will be a cell array of cell arrays of column vectors
dataSeparated = extract_columns(dataOrig, allColNums, 'OutputMode', 'single');

% Initialize a cell array for the processed data
%   Note: each element corresponds to one original column
dataAvgSeparated = cell(size(dataSeparated));

% Extract the column number of interest
vecsOfInterest = dataSeparated(colNumToAverage);

% Average over traces from each group separately
vecsAvg = compute_average_trace(vecsOfInterest, 'Grouping', grouping);

% Store averaged vectors in dataAvgSeparated
dataAvgSeparated(colNumToAverage) = vecsAvg;

% Extract other columns
vecsOther = dataSeparated(colNumOther);

% For other vectors, extract the first vector from each group
vecsFirst = compute_combined_trace(vecsOther, 'first', 'Grouping', grouping);

% Store copied vectors in dataAvgSeparated
dataAvgSeparated(colNumOther) = vecsFirst;

% Concatenate the results from each group
dataAvg = cellfun(@(x) horzcat(x{:}), dataAvgSeparated, ...
                    'UniformOutput', false);

% Reorganize the output so that each element of the cell array are
%   the averaged columns from the same group
dataAvg = extract_columns(dataAvg, allGroupNums, 'OutputMode', 'single');
dataAvg = cellfun(@(x) horzcat(x{:}), dataAvg, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[tVecs, vVecs, iVecs, gVecs] = extract_columns(dataOrig, 1:4);

% For time, current and conductance, 
%   take the first copy and repeat nGroups times
%   Note: Since pulse widths are not necessarily the same, 
%           current vectors should not be averaged
[tVecs, iVecs, gVecs] = ...
    argfun(@(x) repmat(x(1), nGroups, 1), tVecs, iVecs, gVecs);

dataAvg = cellfun(@(x, y, z, w) horzcat(x, y, z, w), ...
                    tVecs, vVecs, iVecs, gVecs, ...
                    'UniformOutput', false);

vecs = dataSeparated{colNumToAverage};

allColNums = create_indices('IndexEnd', nColumns);

% For other vectors, extract the first vector and make nGroups copies
vecsFirst = cellfun(@(x) repmat(x(1), nGroups, 1), vecsOther, ...
                        'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%