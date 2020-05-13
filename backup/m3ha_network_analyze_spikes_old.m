function [oscParams, oscData] = m3ha_network_analyze_spikes (varargin)
%% Analyzes .spi files in a directory
% Usage: [oscParams, oscData] = m3ha_network_analyze_spikes (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       oscParams    - a table with each network as a row and columns:
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
%                   default == fullfile(outFolder, 
%                                   [dirBase, '_oscillation_params.csv'])
%                   - 'PlotFlag': whether to plot analysis results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_files.m
%       cd/apply_over_cells.m
%       cd/argfun.m
%       cd/array_fun.m
%       cd/compute_autocorrelogram.m
%       cd/compute_spike_histogram.m
%       cd/create_subplots.m
%       cd/extract_columns.m
%       cd/extract_fileparts.m
%       cd/load_neuron_outputs.m
%       cd/plot_autocorrelogram.m
%       cd/plot_spike_histogram.m
%       cd/renamevars_custom.m
%       cd/save_all_figtypes.m
%       cd/transpose_table.m
%
% Used by:
%       cd/m3ha_network_launch.m

% File History:
% 2020-01-30 Modified from m3ha_network_plot_essential.m
% 2020-02-05 Added percentActive
% 2020-02-09 Changed definition of hasOscillation
% 2020-02-09 Now sorts by 'datenum'
% 2020-02-09 Add 'PlotFlag' as an optional argument

%% Hard-coded parameters
spiExtension = 'spi';

% File parameters
%   Note: Must be consistent with m3ha_network_launch.m
prefixRT = 'RE';
prefixTC = 'TC';
paramPrefix = 'sim_params';
stimStartStr = 'stimStart';
stimDurStr = 'stimDur';
nCellsStr = 'nCells';
simNumberStr = 'simNumber';

% TODO: Make optional argument
figTypes = 'png';

%% Default values for optional arguments
inFolderDefault = pwd;      % use current directory by default
outFolderDefault = '';      % set later
sheetNameDefault = '';      % no spreadsheet name by default
plotFlagDefault = true;     % plot analysis results by default

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
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
sheetName = iP.Results.SheetName;
plotFlag = iP.Results.PlotFlag;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Preparation
% Set default output folder
if isempty(outFolder)
    outFolder = inFolder;
end

% Look for all RT neuron .spi files
[~, spiPathsRT] = all_files('Directory', inFolder, 'Prefix', prefixRT, ...
                            'Extension', spiExtension, 'SortBy', 'datenum');

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
    sheetName = fullfile(outFolder, [dirBase, '_oscillation_params.csv']);
end

%% Do the job
% Load simulation parameters
paramTables = array_fun(@(x) readtable(x, 'ReadRowNames', true), ...
                        paramFilePaths, 'UniformOutput', false);

% Keep just the Value column in all tables
paramTables = array_fun(@(x) x(:, 'Value'), ...
                        paramTables, 'UniformOutput', false);

% Rename 'Value' by the condition string
paramTables = array_fun(@(x, y) renamevars_custom(x, 'Value', y), ...
                        paramTables, condStr, 'UniformOutput', false);

% Combine all tables
allParamsTable = apply_over_cells(@horzcat, paramTables);

% Initialize the oscillation table with simulation parameters
oscParams = transpose_table(allParamsTable);

% Extract fields
stimStartMs = oscParams.(stimStartStr);
stimDurMs = oscParams.(stimDurStr);
nCells = oscParams.(nCellsStr);
simNumber = oscParams.(simNumberStr);

% Load simulated data
[spikesDataRT, spikesDataTC] = ...
    argfun(@(x) load_neuron_outputs('FileNames', x), spiPathsRT, spiPathsTC);



% Parse spikes
[parsedParams, parsedData] = ...
    cellfun(@(a, b, c, d, e, f) m3ha_network_parse_spikes(a, b, c, d, e, f, ...
                                    plotFlag, outFolder, figTypes), ...
                spikesDataRT, spikesDataTC, num2cell(stimStartMs), ...
                num2cell(stimDurMs), num2cell(nCells), condStr);

% Convert structure arrays to tables
[parsedParamsTable, parsedDataTable] = ...
    argfun(@struct2table, parsedParams, parsedData);

%% Return as output
oscParams = horzcat(oscParams, parsedParamsTable);
oscData = parsedDataTable;

% Reorder according to simNumber
[oscParams, origInd] = sortrows(oscParams, simNumberStr);
oscData = oscData(origInd, :);

%% Save the output
writetable(oscParams, sheetName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
            m3ha_network_parse_spikes (spikesDataRT, spikesDataTC, ...
                                        stimStartMs, stimDurMs, nCells, ...
                                        condStr, plotFlag, outFolder, figTypes)
%% Parse spikes from .spi files

%% Hard-coded parameters
xLimits = stimStartMs/1000 + [0, 10];          % in seconds
figName = '';

% Decide on figure name
if isempty(figName)
    figName = fullfile(outFolder, [condStr, '_analyze_spikes.png']);
end

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

% Put all spike times together
spikeTimesAll = [spikeTimesRT; spikeTimesTC];

% Count the number of active neurons
nActive = numel(unique(cellIdRT)) + numel(unique(cellIdTC));

% Compute the percentage of active neurons
percentActive = 100 * nActive / (nCells * 2);

% Use all spikes to compute an oscillation duration
%   TODO: Modify this for multi-cell layers
[histParams, histData] = ...
    compute_spike_histogram(spikeTimesAll, 'StimStartMs', stimStartMs);

% Use all spikes to compute an oscillation period
%   TODO: Modify this for multi-cell layers
[autoCorrParams, autoCorrData] = ...
    compute_autocorrelogram(spikeTimesAll, 'StimStartMs', stimStartMs, ...
                            'SpikeHistParams', histParams, ...
                            'SpikeHistData', histData);

% Decide whether there is an oscillation based on RT spikes
%   Note: Number of bursts in an oscillation must be more than 2.
hasOscillation = histParams.nBurstsInOsc > 2;

%% Plot for verification
if plotFlag
    % Create a figure
    [fig, ax] = create_subplots(2, 1, 'AlwaysNew', true, ...
                                'FigExpansion', [2, 2]);

    % Add figure title base
    histParams.figTitleBase = replace(condStr, '_', '\_');
    autoCorrParams.figTitleBase = replace(condStr, '_', '\_');

    % Plot spike histogram with burst detection
    subplot(ax(1));
    plot_spike_histogram(histData, histParams, 'XLimits', xLimits);

    % Plot autocorrelation function
    subplot(ax(2));
    plot_autocorrelogram(autoCorrData, autoCorrParams, ...
                            'XLimits', xLimits, ...
                            'PlotType', 'acfFiltered');

    % Save figure
    save_all_figtypes(fig, figName, figTypes);

    % Close figure
    close(fig);
end

%% Save results in output
parsedParams.stimStartMs = stimStartMs;
parsedParams.stimDurMs = stimDurMs;
parsedParams.hasOscillation = hasOscillation;
parsedParams.nActive = nActive;
parsedParams.percentActive = percentActive;
parsedParams = merge_structs(parsedParams, histParams);
parsedParams = merge_structs(parsedParams, autoCorrParams);

parsedData.cellIdRT = cellIdRT;
parsedData.spikeTimesRT = spikeTimesRT;
parsedData.cellIdTC = cellIdTC;
parsedData.spikeTimesTC = spikeTimesTC;
parsedData = merge_structs(parsedData, histData);
parsedData = merge_structs(parsedData, autoCorrData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

allParamsTable = apply_over_cells(@outerjoin, paramTables, 'Keys', 'Row', 'MergeKeys', true);

hasOscillation = ~isempty(spikeTimesTC) && ...
                    any(spikeTimesTC > stimStartMs + stimDurMs);

% Use RT spikes to compute an oscillation duration
[histParams, histData] = ...
    compute_spike_histogram(spikeTimesRT, 'StimStartMs', stimStartMs);
% Use RT spikes to compute an oscillation period
[autoCorrParams, autoCorrData] = ...
    compute_autocorrelogram(spikeTimesRT, 'StimStartMs', stimStartMs, ...
                            'SpikeHistParams', histParams, ...
                            'SpikeHistData', histData);

% Reorder stuff
[simNumber, origInd] = sort(simNumber);
[spikesDataRT, spikesDataTC, stimStartMs, ...
        stimDurMs, nCells, condStr] = ...
    argfun(@(x) x(origInd), ...
            spikesDataRT, spikesDataTC, stimStartMs, ...
            stimDurMs, nCells, condStr);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
