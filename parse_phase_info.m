function varargout = parse_phase_info (pValues, readout, phaseVectors, varargin)
%% Parses phase-related information from parameter and readout values
% Usage: [uniquePhases, phaseBoundaries, averageWindows, phaseAverages, indSelected] = parse_phase_info (pValues, readout, phaseVectors, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       pValues = transpose(1:10);
%       readout = randi(5, 10, 3);
%       phaseVectors = {[1; 1; 2; 2; 2; 2; 2; 2; 3; 3], [1; 1; 1; 1; 2; 2; 3; 3; 3; 3], [1; 1; 1; 1; 1; 1; 2; 2; 3; 3]};
%       TODO: phaseVectors = {[1; 1; 2; 2; 2; 2; 2; 2; 3; 3], [1; 1; 1; 1; 2; 2; 3; 3; 3; 3], [1; 1; 1; 1; 1; 2; 2; 2; 2; 2]};
%       [uniquePhases, phaseBoundaries, averageWindows, phaseAverages, indSelected] = parse_phase_info(pValues, readout, phaseVectors)
%
% Outputs:
%       uniquePhases    - TODO
%                       specified as a TODO
%       phaseBoundaries - TODO
%                       specified as a TODO
%       averageWindows  - TODO
%                       specified as a TODO
%       phaseAverages   - TODO
%                       specified as a TODO
%       indSelected     - TODO
%                       specified as a TODO
%
% Arguments:
%       pValues     - a column vector of parameter values
%                   must be a numeric vector
%       readout     - a readout matrix where each column is a readout vector
%                   must be a numeric matrix or a cell array of numeric vectors
%       phaseVectors- phase vector(s) corresponding to each readout vector
%                   must be a numeric matrix or a cell array of numeric vectors
%       varargin    - 'PhaseBoundaries': vector of phase boundaries
%                   must be a numeric vector
%                   default == set in parse_phase_info.m
%                   - 'AverageWindows': windows to average values
%                       Note: If a matrix cell array, 
%                           each column is for a curve and each row is for a phase
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == set in parse_phase_info.m
%                   - 'PhaseAverages': average values for each phase
%                       Note: If a matrix cell array, 
%                           each column is for a curve and each row is for a phase
%                   must be a numeric 2-D array
%                   default == set in parse_phase_info.m
%                   - 'IndSelected': selected indices to mark differently
%                       Note: If a matrix cell array, 
%                           each column is for a curve and each row is for a phase
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == set in parse_phase_info.m
%                   - 'NLastOfPhase': number of values at the last of a phase
%                   must be empty or a positive integer scalar
%                   default == 10
%                   - 'NToAverage': number of values to average
%                   must be empty or a positive integer scalar
%                   default == set in select_similar_values.m
%                   - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'          - set in select_similar_values.m
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'auto'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be empty or a nonnegative scalar
%                   default == set in select_similar_values.m
%
% Requires:
%       cd/argfun.m
%       cd/compute_phase_average.m
%       cd/compute_index_boundaries.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/force_matrix.m
%       cd/match_format_vector_sets.m
%       cd/unique_custom.m
%       cd/set_default_flag.m
%       cd/trim_nans.m
%
% Used by:
%       cd/plot_bar.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-11-24 Moved from plot_tuning_curve.m
% 2019-11-25 Added phaseBoundaries, averageWindows, etc. as optional arguments
% TODO: Deal with the case when the number of phases are different

%% Hard-coded parameters
validSelectionMethods = {'auto', 'notNaN', 'maxRange2Mean'};
defaultNLastOfPhase = 10;

%% Default values for optional arguments
phaseBoundariesDefault = [];    % set later
averageWindowsDefault = {};     % set later
phaseAveragesDefault = [];      % set later
indSelectedDefault = [];        % set later
nLastOfPhaseDefault = [];       % set later
nToAverageDefault = [];         % set in select_similar_values.m
selectionMethodDefault = 'auto';% set in select_similar_values.m
maxRange2MeanDefault = [];      % set in select_similar_values.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'pValues', ...              % vector of parameter values
    @(x) validateattributes(x, {'numeric', 'datetime', 'duration'}, {'vector'}));
addRequired(iP, 'readout', ...              % a readout matrix
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['readout must be a numeric array ', ...
                    'or a cell array of numeric vectors!']));
addRequired(iP, 'phaseVectors', ...
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['phaseVectors must be a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PhaseBoundaries', phaseBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'AverageWindows', averageWindowsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'PhaseAverages', phaseAveragesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'IndSelected', indSelectedDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'NLastOfPhase', nLastOfPhaseDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['NLastOfPhase must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'NToAverage', nToAverageDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['NToAverage must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) assert(isempty(x) || isscalar(x) && x >= 0, ...
                ['MaxRange2Mean must be either empty ', ...
                    'or a nonnegative scalar!']));

% Read from the Input Parser
parse(iP, pValues, readout, phaseVectors, varargin{:});
phaseBoundaries = iP.Results.PhaseBoundaries;
averageWindows = iP.Results.AverageWindows;
phaseAverages = iP.Results.PhaseAverages;
indSelected = iP.Results.IndSelected;
nLastOfPhase = iP.Results.NLastOfPhase;
nToAverage = iP.Results.NToAverage;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
maxRange2Mean = iP.Results.MaxRange2Mean;

%% Preparation
% Decide on nLastOfPhase
if isempty(nLastOfPhase)
    nLastOfPhase = defaultNLastOfPhase;
end

% Decide whether to compute things
[computePhaseBoundaries, computeAverageWindows, ...
        computePhaseAverages, computeIndSelected] = ...
    argfun(@(x) set_default_flag([], x), ...
            nargout >= 2 && isempty(phaseBoundaries), ...
            nargout >= 3 && isempty(averageWindows), ...
            nargout >= 4 && isempty(phaseAverages), ...
            nargout >= 4 && isempty(indSelected));

% Force as cell arrays of numeric vectors
if ~isempty(phaseBoundaries)
    [readout, phaseBoundaries] = match_format_vector_sets(readout, ...
                                    phaseBoundaries, 'ForceCellOutputs', true);
else
    % Force as cell arrays of numeric vectors
    [readout, phaseVectors] = match_format_vector_sets(readout, ...
                                    phaseVectors, 'ForceCellOutputs', true);
end

% If there is any trailing NaNs, remove them
%   TODO: Make trim_nans accept cell arrays
readout = cellfun(@(x) trim_nans(x, 'trailing'), ...
                    readout, 'UniformOutput', false);

%% Get unique phases
if ~isempty(phaseBoundaries)
    % Count the number of phases for each vector
    nPhases = count_samples(phaseBoundaries) + 1;

    % Set the unique phases to be 1:nPhases
    uniquePhases = arrayfun(@(x) transpose(1:x), ...
                            nPhases, 'UniformOutput', false);
else
    % Get the unique phases across all phase vectors
    uniquePhases = cellfun(@(x) unique_custom(x, 'stable', 'IgnoreNaN', true), ...
                            phaseVectors, 'UniformOutput', false);
end

% Count the number of phases for each readout column
% nPhases = count_samples(uniquePhases);

% Compute the maximum number of phases
% maxNPhases = max(nPhases);

%% Compute phase vectors from boundaries
% if isempty(phaseVectors) && ~isempty(phaseBoundaries)
%     phaseVectors = ...
%         compute_phase_vectors_from_boundaries(pValues, phaseBoundaries);
% end

%% Compute phase boundaries and corresponding indices
if computePhaseBoundaries
    % TODO: Make function [valueBoundaries, indBoundaries] = ...
    %                       compute_value_boundaries(values, grouping)
    % Compute all possible index boundaries for the phases
    indBoundaries = compute_index_boundaries('Grouping', phaseVectors, ...
                                                'TreatNaNsAsGroup', false);

    % Convert to phase units
    phaseBoundaries = extract_subvectors(pValues, 'Indices', indBoundaries);
elseif ~isempty(phaseBoundaries)
    % Compute index boundaries
    indBoundaries = cellfun(@(x) find_closest(pValues, x), phaseBoundaries, ...
                            'UniformOutput', false);
end

%% Compute averaging windows
if computeAverageWindows
    % Compute averaging windows for each set of phase boundaries
    averageWindows = ...
        cellfun(@(x) compute_average_window(pValues, x, nLastOfPhase), ...
                indBoundaries, 'UniformOutput', false);

    % Force as a matrix
    averageWindows = force_matrix(averageWindows, 'TreatCellAsArray', true);
end

%% Compute phase averages and selected indices
% TODO: Use averageWindows if provided
if computePhaseAverages || computeIndSelected        
    % Compute all phase averages for each vector
    %   Note: this generates a cell array of cell arrays of vectors
    [phaseAveragesCellCell, indSelectedCellCell] = ...
        cellfun(@(x, y, z) ...
                compute_all_phase_averages(x, y, z, nLastOfPhase, ...
                                nToAverage, selectionMethod, maxRange2Mean), ...
                readout, indBoundaries, uniquePhases, 'UniformOutput', false);

    % Reorganize so that outputs are a matrix cell array
    %   Note: each row is a phase and each column is a curve
    % TODO: Deal with the case when the number of phases are different
    [phaseAveragesCell, indSelectedCell] = ...
        argfun(@(x) force_matrix(x, 'TreatCellAsArray', true), ...
                phaseAveragesCellCell, indSelectedCellCell);

    % Take scalar phase averages out of the cell array
    phaseAverages = cell2num(phaseAveragesCell);
    indSelected = indSelectedCell;
end

%% Output results
varargout{1} = uniquePhases;
if nargout >= 2
    varargout{2} = phaseBoundaries;
end
if nargout >= 3
    varargout{3} = averageWindows;
end
if nargout >= 4
    varargout{4} = phaseAverages;
end
if nargout >= 5
    varargout{5} = indSelected;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function averageWindow = compute_average_window(pValues, indBoundaries, ...
                                                    nLastOfPhase)
%% Computes windows to compute phase averages

% Compute the last index of each window
iLastEachWindow = [indBoundaries - 0.5; numel(pValues)];

% Compute the first index of each window
iFirstEachWindow = iLastEachWindow - nLastOfPhase + 1;

% Compute each averaging window
averageWindow = ...
    arrayfun(@(u, v) extract_subvectors(pValues, 'Indices', [u; v]), ...
            iFirstEachWindow, iLastEachWindow, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phaseAverages, indSelected] = ...
            compute_all_phase_averages(readout, indBoundaries, uniquePhases, ...
                    nLastOfPhase, nToAverage, selectionMethod, maxRange2Mean)

% Copmute phase average for each unique phase
[phaseAverages, indSelected] = ...
    arrayfun(@(w) compute_phase_average(readout, 'PhaseNumber', w, ...
            'ReturnLastTrial', true, ...
            'PhaseBoundaries', indBoundaries, ...
            'NLastOfPhase', nLastOfPhase, ...
            'NToAverage', nToAverage, ...
            'SelectionMethod', selectionMethod, ...
            'MaxRange2Mean', maxRange2Mean), ...
            uniquePhases, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phaseVectors = ...
                compute_phase_vectors_from_boundaries(pValues, phaseBoundaries)

% TODO
if iscell(phaseBoundaries)
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%