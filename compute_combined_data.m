function [dataAvg, groups] = ...
                compute_combined_data (dataOrig, combineMethod, varargin)
%% Average data according column numbers to average and to a grouping vector
% Usage: [dataAvg, groups] = ...
%               compute_combined_data (dataOrig, combineMethod, varargin)
% Explanation:
%       Given data as a cell array of sweeps, with each
%           sweep's data represented by many columns,
%           this function allows the user to select specific columns
%           to average and to average over specific groups of sweeps.
%       For the columns that are not averaged, the first sweep
%           out of each group is extracted
%       By default, a single set of vectors is considered as one sweep
% Example(s):
%       compute_combined_data({magic(3), magic(3) + 1, magic(3) + 2}, 'mean', 'Grouping', {'a', 'b', 'b'}, 'ColNumToCombine', 2:3)
%       compute_combined_data({magic(3), magic(3) + 1, magic(3) + 2}, 'bootmean', 'Grouping', {'a', 'b', 'b'}, 'ColNumToCombine', 2:3)
% Outputs:
%       dataAvg     - data averaged
%                   specified as a numeric array
%       groups      - unique grouping values
%                   specified as a numeric vector
% Arguments:
%       dataOrig        - original data
%                       must be a numeric array or a cell array of numeric arrays
%       combineMethod   - method for combining traces
%                       see compute_combined_trace.m
%       varargin    - 'ColNumToCombine': column number(s) to average
%                   must be empty or a positive integer vector
%                   default == transpose(1:min(nColumns))
%                   - 'Grouping': a grouping vector used to group traces
%                   must be a vector
%                   default == []
%
% Requires:
%       cd/compute_combined_trace.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/extract_columns.m
%       cd/isnum.m
%       cd/iscellvector.m
%
% Used by:
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-01-12 Adapted from code in m3ha_import_raw_traces.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
colNumToCombineDefault = [];
groupingDefault = [];               % set later

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
addRequired(iP, 'combineMethod');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ColNumToCombine', colNumToCombineDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x), ...
                ['ColNumToCombine must be either empty or , ', ...
                    'a positive integer vector!']));
addParameter(iP, 'Grouping', groupingDefault);

% Read from the Input Parser
parse(iP, dataOrig, combineMethod, varargin{:});
colNumToCombine = iP.Results.ColNumToCombine;
grouping = iP.Results.Grouping;

%% Preparation
% Set default grouping vector if not provided
%   or use it to determine the number of sweeps
if isempty(grouping)
    if isnumeric(dataOrig) || iscellvector(dataOrig)
        % By default, a single set of vectors is considered as one sweep
        nSweeps = 1;
    else
        % Count the number of sweeps
        nSweeps = numel(dataOrig);
    end

    % Make them all the same group
    grouping = ones(nSweeps, 1);
else
    % Count the number of sweeps expected
    nSweeps = numel(grouping);
end

% Count the number of vectors in each cell of dataOrig
nVectors = count_vectors(dataOrig);

% Decide on the number of columns for each sweep
if nSweeps == numel(nVectors)
    % Any single set of vectors is already considered as one sweep
    %   so nVectors is indeed nColumns
    nColumns = nVectors;
    treatCnvAsColumns = true;
elseif nSweeps == nVectors
    % The single set of vectors is considered as multiple sweeps
    %   with one column for each sweep
    nColumns = 1;
    treatCnvAsColumns = false;
end

% Generate a vector of column numbers
allColNums = create_indices('IndexEnd', min(nColumns));

% Set default column numbers to average
if isempty(colNumToCombine)
    % Average all columns by default
    colNumToCombine = allColNums;
end

% Compute other column numbers
colNumOther = setdiff(allColNums, colNumToCombine);

% Find unique grouping values
groups = unique(grouping, 'sorted');

%% Do the job
% Extract each column from all sweeps separately
%   Note: dataSeparated will be a cell array of cell arrays of column vectors
dataSeparated = extract_columns(dataOrig, allColNums, 'OutputMode', 'single', ...
                                'TreatCnvAsColumns', treatCnvAsColumns);

% Initialize a cell array for the processed data
%   Note: each element corresponds to one original column
dataAvgSeparated = cell(size(dataSeparated));

% Extract the column number of interest
vecsOfInterest = dataSeparated(colNumToCombine);

% Average over traces from each group separately
vecsAvg = compute_combined_trace(vecsOfInterest, combineMethod, ...
                                'Grouping', grouping);

% Store averaged vectors in dataAvgSeparated
dataAvgSeparated(colNumToCombine) = vecsAvg;

% Extract other columns
vecsOther = dataSeparated(colNumOther);

% For other vectors, extract the first vector from each group
switch combineMethod
    case {'average', 'mean', 'maximum', 'minimum', ...
            'all', 'any', 'first', 'last'}
        vecsFirst = compute_combined_trace(vecsOther, 'first', ...
                                            'Grouping', grouping);
    case {'bootmean', 'bootmax', 'bootmin'}
        vecsFirst = vecsOther;
    otherwise
        error('combineMethod unrecognized!');
end

% Store copied vectors in dataAvgSeparated
dataAvgSeparated(colNumOther) = vecsFirst;

%  Count the number of vectors in each cell of dataOrig
nOutputs = count_vectors(dataAvgSeparated);

% Generate a vector of output numbers
allOutNums = create_indices('IndexEnd', min(nOutputs));

% Reorganize the output so that each element of the cell array are
%   the averaged columns from the same group
dataAvg = extract_columns(dataAvgSeparated, allOutNums, ...
                                'OutputMode', 'single', ...
                                'TreatCnvAsColumns', true);
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

vecs = dataSeparated{colNumToCombine};

allColNums = create_indices('IndexEnd', nColumns);

% For other vectors, extract the first vector and make nGroups copies
vecsFirst = cellfun(@(x) repmat(x(1), nGroups, 1), vecsOther, ...
                        'UniformOutput', false);

% Concatenate the results from each group
dataAvg = cellfun(@(x) horzcat(x{:}), dataAvgSeparated, ...
                    'UniformOutput', false);
dataAvg = extract_columns(dataAvg, allGroupNums, 'OutputMode', 'single');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
