function [bestParamsTable, bestParamsLabel] = ...
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
%                   - Any other parameter-value pair for 
%                           m3ha_neuron_run_and_analyze()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_label_from_numbers.m
%       cd/extract_fields.m
%       cd/istext.m
%       cd/load_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/set_fields_zero.m
%
% Used by:
%       /media/adamX/m3ha/optimizer4compgabab/singleneuronfitting63.m

% File History:
% 2019-11-23 Created by Adam Lu
% 

%% Hard-coded parameters
validSimModes = {'active', 'passive'};

%% Default values for optional arguments
simModeDefault = 'active';      % simulate active responses by default

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

% Read from the Input Parser
parse(iP, candParamsTablesOrFiles, varargin{:});
simMode = validatestring(iP.Results.SimMode, validSimModes);

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

% Decide on table labels
if isempty(candParamsFiles)
    tableLabels = create_label_from_numbers(1:nTables, 'Prefix', 'table');
else
    tableLabels = candParamsFiles;
end

% Turn off all flags for stats and plots
otherArguments = ...
    set_fields_zero(otherArguments, ...
        'saveLtsInfoFlag', 'saveLtsStatsFlag', ...
        'saveSimCmdsFlag', 'saveStdOutFlag', 'saveSimOutFlag', ...
        'plotConductanceFlag', 'plotCurrentFlag', ...
        'plotIndividualFlag', 'plotResidualsFlag', 'plotOverlappedFlag', ...
        'plotIpeakFlag', 'plotLtsFlag', 'plotStatisticsFlag', ...
        'plotSwpWeightsFlag');

%% Do the job
% Compute errors for all tables
errorStructs = cellfun(@(x) m3ha_neuron_run_and_analyze(x, ...
                        'SimMode', simMode, otherArguments), candParamsTables);

% Extract all total errors
totalErrors = extract_fields(errorStructs, 'totalError', 'UniformOutput', true);

% Find the index of the table with the least error
[totalErrorBest, iTableBest] = min(totalErrors);

%% Output results
% Return the table with the least error
bestParamsTable = candParamsTables{iTableBest};
bestParamsLabel = tableLabels{iTableBest};

% Display result
fprintf('%s has the least error: %g!\n', bestParamsLabel, totalErrorBest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%