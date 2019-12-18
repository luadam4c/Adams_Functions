function [bestParamsTable, bestParamsLabel, errorTable] = ...
                m3ha_neuron_choose_best_params (candParamsTablesOrFiles, varargin)
%% Chooses among candidates the NEURON parameters that fits a cell's data the best
% Usage: [bestParamsTable, bestParamsLabel] = ...
%               m3ha_neuron_choose_best_params (candParamsTablesOrFiles, varargin)
% Explanation:
%       Computes errors for more than one candidate sets of NEURON parameters
%            and choose the one with the least total error as the best 
%
% Example(s):
%       TODO
%
% Outputs:
%       bestParamsTable - the NEURON table for best parameters
%                       specified as a table
%       bestParamsLabel - file name or table name for the best parameters
%                       specified as a character vector
%
% Arguments:
%       candParamsTablesOrFiles  - candidate sets of NEURON parameter
%                                   tables or spreadsheet file names
%                   must be a cell array or string array
%       varargin    - 'SimMode': simulation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - simulate a current pulse response
%                       'active'  - simulate an IPSC response
%                   default == 'active'
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == extract_common_prefix(fileBase)
%                   - Any other parameter-value pair for 
%                           m3ha_neuron_run_and_analyze()
%
% Requires:
%       cd/argfun.m
%       cd/combine_strings.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/extract_fields.m
%       cd/extract_substrings.m
%       cd/isemptycell.m
%       cd/istext.m
%       cd/load_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/set_fields_zero.m
%
% Used by:
%       cd/m3ha_rank_neurons.m
%       /media/adamX/m3ha/optimizer4compgabab/singleneuronfitting63.m

% File History:
% 2019-11-23 Created by Adam Lu
% 2019-11-28 Now saves error table and plots individual plots for each set
%               of parameters
% TODO: Create a combined error file and rank the errors
% 

%% Hard-coded parameters
validSimModes = {'active', 'passive'};
iterStrPattern = 'singleneuronfitting[\d]*';
cellNamePattern = '[A-Z][0-9]{6}';
errorSheetSuffix = 'error_comparison';
errorSheetExtension = 'csv';

%% Default values for optional arguments
simModeDefault = 'active';      % simulate active responses by default
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
prefixDefault = '';             % set later

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
addRequired(iP, 'candParamsTablesOrFiles', ...
    @(x) validateattributes(x, {'cell', 'string'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SimMode', simModeDefault, ...
    @(x) any(validatestring(x, validSimModes)));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, candParamsTablesOrFiles, varargin{:});
simMode = validatestring(iP.Results.SimMode, validSimModes);
outFolder = iP.Results.OutFolder;
prefix = iP.Results.Prefix;

% Keep unmatched arguments for the m3ha_neuron_run_and_analyze() function
otherArguments = iP.Unmatched;

%% Preparation
% Parse first argument
if istext(candParamsTablesOrFiles)
    candParamsFiles = candParamsTablesOrFiles;
    candParamsTables = {};
else
    candParamsTables = candParamsTablesOrFiles;
    candParamsFiles = {};
end

% Load parameters if necessary
if isempty(candParamsTables)
    candParamsTables = cellfun(@load_params, candParamsTablesOrFiles, ...
                                'UniformOutput', false);
end

% Count the number of tables
nTables = numel(candParamsTables);

% If there are no tables, return
if nTables == 0
    bestParamsTable = table.empty;
    bestParamsLabel = '';
    errorTable = table.empty;
    return
end

% Decide on iteration strings and cell names
%   Note: Make sure iterStr and cellName have nTables rows
if isempty(candParamsFiles)
    iterStr = create_labels_from_numbers(1:nTables, 'Prefix', 'table');
    cellName = repmat({'some_cell'}, nTables, 1);
else
    % Extract the chosen iteration string
    iterStr = extract_substrings(candParamsFiles, 'RegExp', iterStrPattern);

    % Extract the cell names
    cellName = extract_substrings(candParamsFiles, 'RegExp', cellNamePattern);
end

% Get unique cell names
uniqueCellNames = unique(cellName);

% Check if all cell names are the same
if numel(uniqueCellNames) > 2
    error('Candidate parameters must all come from the same cell!');
else
    uniqueCellName = uniqueCellNames{1};
end

% Turn off all flags for stats and plots except plotIndividualFlag
% TODO: CHange this to use input parser for case-insensitivity
otherArguments = ...
    set_fields_zero(otherArguments, ...
        'saveLtsInfoFlag', 'saveLtsStatsFlag', ...
        'saveSimCmdsFlag', 'saveStdOutFlag', 'saveSimOutFlag', ...
        'plotConductanceFlag', 'plotCurrentFlag', ...
        'plotResidualsFlag', 'plotOverlappedFlag', ...
        'plotIpeakFlag', 'plotLtsFlag', 'plotStatisticsFlag', ...
        'plotSwpWeightsFlag');
otherArguments = ...
    set_fields_zero(otherArguments, ...
        'SaveLtsInfoFlag', 'SaveLtsStatsFlag', ...
        'SaveSimCmdsFlag', 'SaveStdOutFlag', 'SaveSimOutFlag', ...
        'PlotConductanceFlag', 'PlotCurrentFlag', ...
        'PlotResidualsFlag', 'PlotOverlappedFlag', ...
        'PlotIpeakFlag', 'PlotLtsFlag', 'PlotStatisticsFlag', ...
        'PlotSwpWeightsFlag');

% Decide on prefix if not provided
if isempty(prefix)
    % Use the cell name
    prefix = uniqueCellName;

    % If no cell name, use the directory base
    if isempty(prefix)
        prefix = extract_fileparts(outFolder, 'dirbase');
    end
end

% Create candidate labels
candLabel = combine_strings('Substrings', {prefix, 'from', iterStr});

%% Do the job
% Compute errors for all tables
errorStructs = cellfun(@(x, y) m3ha_neuron_run_and_analyze(x, ...
                            'PlotIndividualFlag', true, ...
                            'SimMode', simMode, 'OutFolder', outFolder, ...
                            'Prefix', y, otherArguments), ...
                            candParamsTables, candLabel);

% Extract scalar fields of interest
%   Note: must be consistent with compute_single_neuron_errors.m
[totalError, lts2SweepErrorRatio, avgSwpError, avgLtsError, ...
        ltsMatchError, missedLtsError, falseLtsError, ...
        avgLtsAmpError, avgLtsDelayError, avgLtsSlopeError] = ...
    argfun(@(x) extract_fields(errorStructs, x, 'UniformOutput', true), ...
            'totalError', 'lts2SweepErrorRatio', 'avgSwpError', ...
            'avgLtsError', 'ltsMatchError', ...
            'missedLtsError', 'falseLtsError', ...
            'avgLtsAmpError', 'avgLtsDelayError', 'avgLtsSlopeError');

% Extract vector fields of interest
%   Note: must be consistent with compute_single_neuron_errors.m
[ltsFeatureWeights, swpErrors, ltsAmpErrors, ...
        ltsDelayErrors, ltsSlopeErrors] = ...
    argfun(@(x) extract_fields(errorStructs, x, 'UniformOutput', false), ...
            'ltsFeatureWeights', 'swpErrors', 'ltsAmpErrors', ...
            'ltsDelayErrors', 'ltsSlopeErrors');

% Find the index of the table with the least error
[totalErrorBest, iTableBest] = min(totalError);

% Add variables in the beginning
errorTable = table(candLabel, cellName, iterStr, ...
                    totalError, lts2SweepErrorRatio, ...
                    avgSwpError, avgLtsError, ltsMatchError, ...
                    missedLtsError, falseLtsError, ltsFeatureWeights, ...
                    avgLtsAmpError, avgLtsDelayError, avgLtsSlopeError, ...
                    swpErrors, ltsAmpErrors, ltsDelayErrors, ltsSlopeErrors, ...
                    'RowNames', iterStr);

%% Save results
% Create full path to error sheet file
sheetName = strcat(prefix, '_', errorSheetSuffix, '.', errorSheetExtension);
sheetPath = fullfile(outFolder, sheetName);

% Save the error table
writetable(errorTable, sheetPath);

%% Output results
% Return the table with the least error
bestParamsTable = candParamsTables{iTableBest};
bestParamsLabel = candLabel{iTableBest};

% Display result
fprintf('%s has the least error: %g!\n', bestParamsLabel, totalErrorBest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Convert error struct array to a table
errorTable = struct2table(errorStructs, 'AsArray', true);
% Make iterStr row names
errorTable.Properties.RowNames = iterStr;
% Add variables in the beginning
errorTable = addvars(errorTable, candLabel, cellName, ...
                        iterStr, 'Before', 1);
% Create candidate labels
candLabel = strcat(cellName, '_from_', iterStr);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
