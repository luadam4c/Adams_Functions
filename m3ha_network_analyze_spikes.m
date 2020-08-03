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
%       oscParams    - a table with each network as a row 
%                       (using condStr as a row name) and columns:
%                       stimStartMs
%                       stimDurMs
%                       hasOscillation
%                       nActive
%                       nActiveRT
%                       nActiveTC
%                       percentActive
%                       percentActiveRT
%                       percentActiveTC
%                       halfActiveTimeMsRT
%                       halfActiveTimeMsTC
%                       halfActiveLatencyMsRT
%                       halfActiveLatencyMsTC
%                       maxOpenProbabilityDiscrepancy
%                       maxLogOpenProbabilityDiscrepancy
%                       passedOpdThreshold
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
%                   - 'MatFileName': matlab file path for saving results
%                   must be a string scalar or a character vector
%                   default == fullfile(outFolder, 
%                                   [dirBase, '_oscillation_data.mat'])
%                   - 'PlotFlag': whether to plot analysis results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/combine_param_tables.m
%       cd/compute_autocorrelogram.m
%       cd/compute_spike_histogram.m
%       cd/create_subplots.m
%       cd/extract_columns.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
%       cd/read_neuron_outputs.m
%       cd/plot_autocorrelogram.m
%       cd/plot_spike_histogram.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/m3ha_network_launch.m
%       cd/m3ha_plot_figure07.m

% File History:
% 2020-01-30 Modified from m3ha_network_plot_essential.m
% 2020-02-05 Added percentActive
% 2020-02-09 Changed definition of hasOscillation
% 2020-02-09 Now sorts by 'datenum'
% 2020-02-09 Add 'PlotFlag' as an optional argument
% 2020-04-08 Now removes spike times less than stimulation start
% 2020-04-08 Added 'MatFileName' as an optional argument
% 2020-04-08 Now computes and plots precent activated over time
% 2020-04-13 Now computes half rise time
% 2020-05-18 Changed binWidthMs from 10 ms to 100 ms
% 2020-05-18 Added minRelSecProm and made it 0.5
% 2020-08-01 Now analyzes maximum open probability discrepancy
%               if .singsp files are present
% 2020-08-02 Added passedOpdThreshold
%               

%% Hard-coded parameters
spiExtension = 'spi';
specialExtension = 'singsp';

% For compute_spike_histogram.m
binWidthMs = 100;               % use a bin width of 100 ms for simulated networks
minBurstLengthMs = 60;          % bursts must be at least 60 ms
maxFirstInterBurstIntervalMs = 2000;
maxInterBurstIntervalMs = 2000; % bursts are no more than 
                                %   2 seconds apart
minSpikeRateInBurstHz = 100;    % bursts must have a spike rate of 
                                %   at least 100 Hz by default

% For compute_autocorrelogram.m
filterWidthMs = 100;
minRelProm = 0.02;
minRelSecProm = 0.5;            % for simulated networks only

% For compute_activation_profile.m
actProfileBinWidthMs = 500;

% File parameters
%   Note: Must be consistent with m3ha_network_launch.m
prefixRT = 'RE';
prefixTC = 'TC';
prefixTC0 = 'TC[0]';
paramPrefix = 'sim_params';
stimStartStr = 'stimStart';
stimDurStr = 'stimDur';
tStopStr = 'tStop';
nCellsStr = 'nCells';
simNumberStr = 'simNumber';

% TODO: Make optional argument
verbose = true;
figTypes = 'png';

%% Default values for optional arguments
inFolderDefault = pwd;      % use current directory by default
outFolderDefault = '';      % set later
sheetNameDefault = '';      % no spreadsheet name by default
matFileNameDefault = '';    % no mat file name by default
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
addParameter(iP, 'MatFileName', matFileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
sheetName = iP.Results.SheetName;
matFileName = iP.Results.MatFileName;
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

% If nothing found, return
if isempty(spiPathsRT)
    oscParams = table.empty;
    oscData = table.empty;
    return
end

% Create paths to corresponding params files
spiPathBasesRT = extractBefore(spiPathsRT, ['.', spiExtension]);
paramFileBases = replace(spiPathBasesRT, prefixRT, paramPrefix);
paramFilePaths = strcat(paramFileBases, '.csv');

% Create paths to corresponding TC spi files
spiFileBasesTC = replace(spiPathBasesRT, prefixRT, prefixTC);
spiPathsTC = strcat(spiFileBasesTC, ['.', spiExtension]);

% Extract the condition strings
condStr = extractAfter(spiPathBasesRT, [prefixRT, '_']);

% Look for the TC[0] special file corresponding to this condition string
[~, specialPathsTC] = ...
    find_matching_files(condStr, 'ReturnEmpty', true, 'Directory', inFolder, ...
                        'Prefix', prefixTC0, 'Extension', specialExtension);

% Extract the directory base
dirBase = extract_fileparts(inFolder, 'dirbase');

% Decide on spreadsheet name
if isempty(sheetName)
    sheetName = fullfile(outFolder, [dirBase, '_oscillation_params.csv']);
end

% Decide on matfile name
if isempty(matFileName)
    matFileName = fullfile(outFolder, [dirBase, '_oscillation_data.mat']);
end

%% Do the job
% Print message
if verbose
    fprintf('Analyzing Spikes for %s ... \n', dirBase);
end

% Initialize the oscillation table with simulation parameters
oscParams = combine_param_tables(paramFilePaths, 'NewRowNames', condStr);

% Extract fields
stimStartMs = oscParams.(stimStartStr);
stimDurMs = oscParams.(stimDurStr);
tStopMs = oscParams.(tStopStr);
nCells = oscParams.(nCellsStr);
simNumber = oscParams.(simNumberStr);

% Load simulated data
[spikesDataRT, spikesDataTC] = ...
    argfun(@(x) read_neuron_outputs('FileNames', x), spiPathsRT, spiPathsTC);

% Parse spikes
[parsedParams, parsedData] = ...
    cellfun(@(a, b, c, d, e, f, g, h) ...
                    m3ha_network_parse_spikes(a, b, c, d, e, f, g, h, ...
                            plotFlag, outFolder, figTypes, ...
                            binWidthMs, minBurstLengthMs, ...
                            maxFirstInterBurstIntervalMs, ...
                            maxInterBurstIntervalMs, minSpikeRateInBurstHz, ...
                            filterWidthMs, minRelProm, minRelSecProm, ...
                            actProfileBinWidthMs), ...
                spikesDataRT, spikesDataTC, ...
                num2cell(stimStartMs), num2cell(stimDurMs), ...
                num2cell(tStopMs), num2cell(nCells), condStr, specialPathsTC);

% Convert structure arrays to tables
[parsedParamsTable, parsedDataTable] = ...
    argfun(@struct2table, parsedParams, parsedData);

%% Return as output
oscParams = horzcat(oscParams, parsedParamsTable);
oscData = parsedDataTable;
oscData.Properties.RowNames = condStr;

% Reorder according to simNumber
[oscParams, origInd] = sortrows(oscParams, simNumberStr);
oscData = oscData(origInd, :);

%% Save the output
writetable(oscParams, sheetName, 'WriteRowNames', true);
save(matFileName, 'oscParams', 'oscData', '-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
            m3ha_network_parse_spikes (spikesDataRT, spikesDataTC, ...
                stimStartMs, stimDurMs, tStopMs, nCells, ...
                condStr, specialPathTC, plotFlag, outFolder, figTypes, ...
                binWidthMs, minBurstLengthMs, maxFirstInterBurstIntervalMs, ...
                maxInterBurstIntervalMs, minSpikeRateInBurstHz, ...
                filterWidthMs, minRelProm, minRelSecProm, actProfileBinWidthMs)
%% Parse spikes from .spi files

%% Hard-coded parameters
MS_PER_S = 1000;
xLimitsAcf = [0, 10];          % in seconds
xLimitsHist = stimStartMs/1000 + xLimitsAcf;          % in seconds
figName = '';
itm2hDiffLowerLimit = 1e-9;

% Column numbers for simulated data
%   Note: Must be consistent with m3ha_net.hoc
TC_TIME = 1;
TC_VOLT = 2;
TC_IT_M_DEND2 = 10;
TC_IT_MINF_DEND2 = 11;
TC_IT_H_DEND2 = 12;
TC_IT_HINF_DEND2 = 13;

% Compute stimulation end
stimEndMs = stimStartMs + stimDurMs;

% Compute time bins for activation profiles
nBins = floor(tStopMs/actProfileBinWidthMs);
actProfileTimeBinsMs = create_time_vectors(nBins, 'TimeUnits', 'ms', ...
                            'SamplingIntervalMs', actProfileBinWidthMs);

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
[cellIdRT, spikeTimesMsRT] = ...
    extract_columns(spikesDataRT, [RT_CELLID, RT_SPIKETIME]);
[cellIdTC, spikeTimesMsTC] = ...
    extract_columns(spikesDataTC, [TC_CELLID, TC_SPIKETIME]);

% Remove spike times less than stimulation end
toRemoveRT = spikeTimesMsRT < stimEndMs;
cellIdRT(toRemoveRT) = [];
spikeTimesMsRT(toRemoveRT) = [];
toRemoveTC = spikeTimesMsTC < stimEndMs;
cellIdTC(toRemoveTC) = [];
spikeTimesMsTC(toRemoveTC) = [];

% Put all spike times together
spikeTimesAll = [spikeTimesMsRT; spikeTimesMsTC];

% Count the number of active neurons
nActiveRT = numel(unique(cellIdRT));
nActiveTC = numel(unique(cellIdTC));
nActive = nActiveRT + nActiveTC;

% Compute the percentage of active neurons
percentActiveRT = 100 * nActiveRT ./ nCells;
percentActiveTC = 100 * nActiveTC ./ nCells;
percentActive = 100 * nActive ./ (nCells * 2);

% Compute the percentage of active neurons over time
percentActivatedRT = compute_activation_profile(cellIdRT, spikeTimesMsRT, ...
                            'TimeBins', actProfileTimeBinsMs, 'NCells', nCells);
percentActivatedTC = compute_activation_profile(cellIdTC, spikeTimesMsTC, ...
                            'TimeBins', actProfileTimeBinsMs, 'NCells', nCells);

% Compute the time of half activation
%   Note: This must be after stimEndMs
halfActiveTimeMsRT = compute_half_rise_time(actProfileTimeBinsMs, ...
                                        percentActivatedRT, percentActiveRT);
halfActiveTimeMsTC = compute_half_rise_time(actProfileTimeBinsMs, ...
                                        percentActivatedTC, percentActiveTC);
halfActiveTimeMsRT(halfActiveTimeMsRT < stimEndMs) = stimEndMs;
halfActiveTimeMsTC(halfActiveTimeMsTC < stimEndMs) = stimEndMs;

% Compute the latency to half activation
halfActiveLatencyMsRT = halfActiveTimeMsRT - stimEndMs;
halfActiveLatencyMsTC = halfActiveTimeMsTC - stimEndMs;

% Store time bins in seconds
[timeBinsSeconds, stimEndSeconds, ...
        halfActiveTimeSecondsRT, halfActiveTimeSecondsTC] = ...
    argfun(@(x) x ./ MS_PER_S, ...
            actProfileTimeBinsMs, stimEndMs, ...
            halfActiveTimeMsRT, halfActiveTimeMsTC);

% Use all spikes to compute an oscillation duration
[histParams, histData] = ...
    compute_spike_histogram(spikeTimesAll, 'StimStartMs', stimStartMs, ...
            'BinWidthMs', binWidthMs, 'MinBurstLengthMs', minBurstLengthMs, ...
            'MaxFirstInterBurstIntervalMs', maxFirstInterBurstIntervalMs, ...
            'MaxInterBurstIntervalMs', maxInterBurstIntervalMs, ...
            'MinSpikeRateInBurstHz', minSpikeRateInBurstHz);

% Use all spikes to compute an oscillation period
[autoCorrParams, autoCorrData] = ...
    compute_autocorrelogram(spikeTimesAll, 'StimStartMs', stimStartMs, ...
                            'SpikeHistParams', histParams, ...
                            'SpikeHistData', histData, ...
                            'FilterWidthMs', filterWidthMs, ...
                            'MinRelProm', minRelProm, ...
                            'MinRelSecProm', minRelSecProm);

% Decide whether there is an oscillation based on all spikes
%   Note: Number of bursts in an oscillation must be more than 2.
%           i.e., there must be at least 2 cycles
hasOscillation = histParams.nBurstsInOsc > 2;

% If a special TC neuron file exists 
if ~isempty(specialPathTC)
    % Read simulation outputs starting from stimulation start
    simDataTC = read_neuron_outputs('FileNames', specialPathTC, ...
                                    'TimeWindows', [stimStartMs; Inf]);

    % Extract columns
    [itmDend2, itminfDend2, ithDend2, ithinfDend2] = ...
        extract_columns(simDataTC, [TC_IT_M_DEND2, TC_IT_MINF_DEND2, ...
                                    TC_IT_H_DEND2, TC_IT_HINF_DEND2]);

    % Compute m2h and its steady state
    itm2h = (itmDend2 .^ 2) .* ithDend2;
    itminf2hinf = (itminfDend2 .^ 2) .* ithinfDend2;

    % Compute m2h discrepancy
    itm2hDiff = itm2h - itminf2hinf;
    itm2hDiff(itm2hDiff < itm2hDiffLowerLimit) = itm2hDiffLowerLimit;

    % Compute maximum m2h discrepancy
    maxOpenProbabilityDiscrepancy = max(itm2hDiff);
    maxLogOpenProbabilityDiscrepancy = log10(maxOpenProbabilityDiscrepancy);
    passedOpdThreshold = maxLogOpenProbabilityDiscrepancy >= -2;
else
    maxOpenProbabilityDiscrepancy = NaN;
    maxLogOpenProbabilityDiscrepancy = NaN;
    passedOpdThreshold = NaN;
end

%% Plot for verification
if plotFlag
    % Create a figure
    [fig, ax] = create_subplots(3, 1, 'AlwaysNew', true, ...
                                'FigExpansion', [2, 3]);

    % Add figure title base
    figTitleBase = replace(condStr, '_', '\_');
    histParams.figTitleBase = figTitleBase;
    autoCorrParams.figTitleBase = figTitleBase;

    % Plot activation profile
    subplot(ax(1)); hold on
    p1 = plot(timeBinsSeconds, percentActivatedRT, 'r', ...
            'DisplayName', 'RT', 'LineWidth', 1);
    p2 = plot(timeBinsSeconds, percentActivatedTC, 'g', ...
            'DisplayName', 'TC', 'LineWidth', 1);
    ylim([0, nCells]);
    v1 = plot_vertical_line(halfActiveTimeSecondsRT, ...
                    'ColorMap', 'r', 'LineStyle', '--');
    v2 = plot_vertical_line(halfActiveTimeSecondsTC, ...
                    'ColorMap', 'g', 'LineStyle', '--');
    v3 = plot_vertical_line(stimEndSeconds, ...
                    'ColorMap', 'k', 'LineStyle', '--');
    xlabel('Time (s)');
    ylabel('Percent Activated (%)');
    title(['Activation profile for ', figTitleBase]);
    legend([p1, p2], 'location', 'northeast');

    % Plot spike histogram with burst detection
    subplot(ax(2));
    plot_spike_histogram(histData, histParams, 'XLimits', xLimitsHist);

    % Plot autocorrelation function
    subplot(ax(3));
    plot_autocorrelogram(autoCorrData, autoCorrParams, ...
                            'XLimits', xLimitsAcf, 'PlotType', 'acfFiltered');

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
parsedParams.nActiveRT = nActiveRT;
parsedParams.nActiveTC = nActiveTC;
parsedParams.percentActive = percentActive;
parsedParams.percentActiveRT = percentActiveRT;
parsedParams.percentActiveTC = percentActiveTC;
parsedParams.halfActiveTimeMsRT = halfActiveTimeMsRT;
parsedParams.halfActiveTimeMsTC = halfActiveTimeMsTC;
parsedParams.halfActiveLatencyMsRT = halfActiveLatencyMsRT;
parsedParams.halfActiveLatencyMsTC = halfActiveLatencyMsTC;
parsedParams.maxOpenProbabilityDiscrepancy = maxOpenProbabilityDiscrepancy;
parsedParams.maxLogOpenProbabilityDiscrepancy = maxLogOpenProbabilityDiscrepancy;
parsedParams.passedOpdThreshold = passedOpdThreshold;
parsedParams = merge_structs(parsedParams, histParams);
parsedParams = merge_structs(parsedParams, autoCorrParams);

parsedData.cellIdRT = cellIdRT;
parsedData.spikeTimesMsRT = spikeTimesMsRT;
parsedData.cellIdTC = cellIdTC;
parsedData.spikeTimesMsTC = spikeTimesMsTC;
parsedData.timeBinsSeconds = timeBinsSeconds;
parsedData.percentActivatedRT = percentActivatedRT;
parsedData.percentActivatedTC = percentActivatedTC;
parsedData = merge_structs(parsedData, histData);
parsedData = merge_structs(parsedData, autoCorrData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function halfRiseTime = compute_half_rise_time (tVec, yVec, yMax)
% TODO: Pull out as its own function

halfRiseTimeRight = tVec(find(yVec >= yMax / 2, 1));
halfRiseTimeLeft = tVec(find(yVec < yMax / 2, 1, 'last'));
halfRiseTime = mean([halfRiseTimeRight, halfRiseTimeLeft]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

allParamsTable = apply_over_cells(@outerjoin, paramTables, 'Keys', 'Row', 'MergeKeys', true);

hasOscillation = ~isempty(spikeTimesMsTC) && ...
                    any(spikeTimesMsTC > stimStartMs + stimDurMs);

% Use RT spikes to compute an oscillation duration
[histParams, histData] = ...
    compute_spike_histogram(spikeTimesMsRT, 'StimStartMs', stimStartMs);
% Use RT spikes to compute an oscillation period
[autoCorrParams, autoCorrData] = ...
    compute_autocorrelogram(spikeTimesMsRT, 'StimStartMs', stimStartMs, ...
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
