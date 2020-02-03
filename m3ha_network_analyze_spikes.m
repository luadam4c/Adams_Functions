function oscTable = m3ha_network_analyze_spikes (varargin)
%% Analyzes .spi files in a directory
% Usage: oscTable = m3ha_network_analyze_spikes (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       oscTable    - a table with each network as a row and columns:
%                       condStr     - condition string
%                   specified as a 2D table
%
% Arguments:
%       varargin    - 'InFolder': directory containing the .singsp files
%                   must be a string scalar or a character vector
%                   default == pwds
%                   - 'OutFolder': output folder
%                   must be a string scalar or a character vector
%                   default == inFolder
%                   - 'SheetName': spreadsheet path for saving
%                   must be a string scalar or a character vector
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_files.m
%       cd/apply_over_cells.m
%       cd/argfun.m
%       cd/array_fun.m
%       cd/extract_columns.m
%       cd/extract_fileparts.m
%       cd/load_neuron_outputs.m
%       cd/renamevars.m
%       cd/transpose_table.m
%
% Used by:

% File History:
% 2020-01-30 Modified from m3ha_network_plot_essential.m

%% Hard-coded parameters
spiExtension = 'spi';

% File parameters
%   Note: Must be consistent with m3ha_network_launch.m
prefixRT = 'RE';
prefixTC = 'TC';
paramPrefix = 'sim_params';
stimStartStr = 'stimStart';
stimDurStr = 'stimDur';

% TODO: Make optional arguments
figTypes = 'png';

%% Default values for optional arguments
inFolderDefault = pwd;      % use current directory by default
outFolderDefault = '';      % set later
sheetNameDefault = '';      % no spreadsheet name by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
sheetName = iP.Results.SheetName;

% Keep unmatched arguments for the TODO() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default output folder
if isempty(outFolder)
    outFolder = inFolder;
end

% Look for all RT neuron .spi files
[~, spiPathsRT] = all_files('Directory', inFolder, 'Prefix', prefixRT, ...
                            'Extension', spiExtension);

% Create paths to corresponding params files
spiPathBasesRT = extractBefore(spiPathsRT, ['.', spiExtension]);
paramFileBases = replace(spiPathBasesRT, prefixRT, paramPrefix);
paramFilePaths = strcat(paramFileBases, '.csv');

% Create paths to corresponding TC spi files
spiFileBasesTC = replace(spiPathBasesRT, prefixRT, prefixTC);
spiPathsTC = strcat(spiFileBasesTC, ['.', spiExtension]);

% Extract the condition strings
condStr = extractAfter(spiPathBasesRT, [prefixRT, '_']);

% Decide on spreadsheet name
if isempty(sheetName)
    dirBase = extract_fileparts(inFolder, 'dirbase');
    sheetName = fullfile(outFolder, [dirBase, '_oscillation_properties.csv']);
end

%% Do the job
% Load simulation parameters
paramTables = array_fun(@(x) readtable(x, 'ReadRowNames', true), ...
                        paramFilePaths, 'UniformOutput', false);

% Keep just the Value column in all tables
paramTables = array_fun(@(x) x(:, 'Value'), ...
                        paramTables, 'UniformOutput', false);

% Rename 'Value' by the condition string
paramTables = array_fun(@(x, y) renamevars(x, 'Value', y), ...
                        paramTables, condStr, 'UniformOutput', false);

% Combine all tables
allParamsTable = apply_over_cells(@outerjoin, paramTables, 'Keys', 'Row', 'MergeKeys', true);

% Initialize the oscillation table with simulation parameters
oscTable = transpose_table(allParamsTable);

% Extract the stimulation start time
stimStartMs = oscTable.(stimStartStr);
stimDurMs = oscTable.(stimDurStr);

% Load simulated data
[spikesDataRT, spikesDataTC] = ...
    argfun(@(x) load_neuron_outputs('FileNames', x), spiPathsRT, spiPathsTC);

% Parse spikes
[parsedParams, parsedData] = ...
    cellfun(@(a, b, c, d) m3ha_network_parse_spikes(a, b, c, d), ...
                spikesDataRT, spikesDataTC, num2cell(stimStartMs), ...
                num2cell(stimDurMs));

%% Add variable to the table
% oscTable = addvars(oscTable, )

%% Save the output
writetable(oscTable, sheetName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
            m3ha_network_parse_spikes (spikesDataRT, spikesDataTC, ...
                                        stimStartMs, stimDurMs)
%% Parse spikes from .spi files

% Column numbers for .spi files
%   Note: Must be consistent with m3ha_net.hoc
RT_CELLID = 1;
RT_SPIKETIME = 2;

TC_CELLID = 1;
TC_SPIKETIME = 2;

% Extract vectors from simulated data
[cellIdRT, spikeTimesRT] = ...
    extract_columns(spikesDataRT, [RT_CELLID, RT_SPIKETIME]);
[cellIdTC, spikeTimesTC] = ...
    extract_columns(spikesDataTC, [TC_CELLID, TC_SPIKETIME]);

% Determine the start of detection
detectStartMs = stimStartMs + stimDurMs;

% Decide whether there is an oscillation
hasOscillation = ~isempty(spikeTimesTC) && any(spikeTimesTC > detectStartMs);

% Compute an oscillation duration
% TODO

% Compute an oscillation period
% TODO

%% Save results in output
parsedParams.stimStartMs = stimStartMs;
parsedParams.stimDurMs = stimDurMs;
parsedParams.oscDurationSec = oscDurationSec;
parsedParams.oscPeriodMs = oscPeriodMs;

parsedData.cellIdRT = cellIdRT;
parsedData.spikeTimesRT = spikeTimesRT;
parsedData.cellIdTC = cellIdTC;
parsedData.spikeTimesTC = spikeTimesTC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
