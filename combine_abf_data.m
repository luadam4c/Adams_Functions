function allData = combine_abf_data (abfPaths, varargin)
%% Combine data from many .abf files and return a structure
% Usage: allData = combine_abf_data (abfPaths, varargin)
% Explanation:
%       TODO
%       Note: Current and voltage vectors are identified using 
%               identify_channels.m by default. If it's already labelled
%               correctly in the abf files, set 'UseOriginal' to be true.
%
% Example(s):
%       [~, abfPaths] = all_files('Ext', 'abf', 'SortBy', 'date');
%       allData = combine_abf_data(abfPaths, 'SaveMatFlag', false)
%       combine_abf_data(abfPaths);
%       combine_abf_data(abfPaths, 'UseOriginal', true);
%       combine_abf_data(abfPaths, 'SaveSheetFlag', false);
%
% Outputs:
%       allData     - a structure with fields:
%                       siMs        - 
%                       vVecs       - voltage vectors
%                       iVecs       - current vectors
%                       phaseBoundaries - 
%                       phaseStrs   - 
%                   specified as a scalar structure
%
% Arguments:
%       abfPaths    - full paths to abf files to combine
%                   must be a character vector, a string array 
%                       or a cell array of character arrays
%       varargin    - 'SaveMatFlag': whether to save combined data as mat file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'RegexpPhaseStr': phase string regular expression
%                   must be a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for parse_all_abfs() 
%                       or parse_abf()
%
% Requires:
%       cd/argfun.m
%       cd/combine_data_from_same_slice.m
%       cd/compute_index_boundaries.m
%       cd/count_A_each_C.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/extract_fileparts.m
%       cd/force_matrix.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/combine_data_from_same_slice.m

% File History:
% 2019-08-23 Pulled from combine_data_from_same_slice.m
% 2019-08-23 Now passes unmatched optional arguments to parse_all_abfs
% TODO: Combine gVecs as well
% 

%% Hard-coded parameters

%% Default values for optional arguments
saveMatFlagDefault = true;      % save combined data by default
regexpPhaseStrDefault = '';     % no phase string regular expression by default

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
addRequired(iP, 'abfPaths', ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RegexpPhaseStr', regexpPhaseStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, abfPaths, varargin{:});
saveMatFlag = iP.Results.SaveMatFlag;
regexpPhaseStr = iP.Results.RegexpPhaseStr;

% Keep unmatched arguments for the parse_all_abfs() or parse_abf() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Extract file bases
allFileNames = extract_fileparts(abfPaths, 'name');

% Extract phase strings
if isempty(regexpPhaseStr)
    allPhaseStrs = extract_fileparts(allFileNames, 'distinct', ...
                                    'Delimiter', '_');
else
    % Note: must pass in file names with the extension
    %           or else extract_fileparts() thinks they are directories
    allPhaseStrs = extract_fileparts(allFileNames, 'base', ...
                                    'RegExp', regexpPhaseStr);
end

% Get the unique phase strings in original order
phaseStrs = unique(allPhaseStrs, 'stable');

% Make sure it is a cell array
phaseStrs = force_column_cell(phaseStrs);

% Count the number of phases
nPhases = numel(phaseStrs);

%% Extract data to combine
% Parse all multi-unit recordings for this slice
[allParams, allAbfData] = ...
    parse_all_abfs('FileNames', abfPaths, otherArguments{:});

% Extract parameters, then clear unused parameters
siMs = allParams.siMs;
clear allParams;

% Extract data, then clear unused data
vVecsAll = allAbfData.vVecs;
iVecsAll = allAbfData.iVecs;
clear allAbfData;

% Find the indices for each phase
indEachPhase = cellfun(@(x) find_in_strings(x, allFileNames), ...
                        phaseStrs, 'UniformOutput', false);

%% Order the data correctedly (may not be needed)
% Put them all together
sortOrder = vertcat(indEachPhase{:});

% Reorder data so that the order matches that of phaseStrs
[siMsSorted, vVecsSorted, iVecsSorted] = ...
    argfun(@(x) x(sortOrder), siMs, vVecsAll, iVecsAll);

%% Combine the data
% Compute the new siMs
siMsSl = mean(siMsSorted);

% Concatenate vectors
% TODO: Fix force_matrix.m to accept cell arrays of non-vectors with an
% optional argument
% [vVecsCombined, iVecsCombined] = argfun(@force_matrix, vVecsSorted, iVecsSorted);
[vVecsCombined, iVecsCombined] = argfun(@(x) horzcat(x{:}), vVecsSorted, iVecsSorted);

%% Create phase boundaries
% Count the number of phase boundaries
nBoundaries = nPhases - 1;

% Compute phase boundaries
if nBoundaries == 0
    phaseBoundaries = [];
else
    % Count the number of sweeps in each file
    nSweepsEachFile = cellfun(@count_vectors, vVecsSorted);

    % Count the number of files for each phase
    nFilesEachPhase = cellfun(@count_samples, indEachPhase);

    % Count the number of sweeps in each phase
    nSweepsEachPhase = count_A_each_C(nSweepsEachFile, nFilesEachPhase);

    % Compute the phase boundaries
    phaseBoundaries = compute_index_boundaries('NEachGroup', nSweepsEachPhase);
end

%% Output as a structure
allData.vVecs = vVecsCombined;
allData.siMs = siMsSl;
allData.iVecs = iVecsCombined;
allData.phaseBoundaries = phaseBoundaries;
allData.phaseStrs = phaseStrs;

%% Save to a matfile if requested
if saveMatFlag
    % Create a matfile path
    commonPrefix = extract_fileparts(abfPaths, 'commonprefix');
    commonDir = extract_fileparts(abfPaths, 'commondirectory');
    matPath = fullfile(commonDir, [commonPrefix, '.mat']);

    % Save data for this slice
    save(matPath, '-struct', 'allData');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%