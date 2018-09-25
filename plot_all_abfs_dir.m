function [abfParamsAllStruct, dataAll, tVecAll, vVecsAll, iVecsAll, ...
            gVecsAll, dataReorderedAll, abfParamsAllCell] = ...
                plot_all_abfs_dir (varargin)
%% Plots all abf files in a directory
% Usage: [abfParamsAllStruct, dataAll, tVecAll, vVecsAll, iVecsAll, ...
%           gVecsAll, dataReorderedAll, abfParamsAllCell] = ...
%               plot_all_abfs_dir (varargin)
%
% Arguments:
%       varargin    - 'Directory': the name of the directory containing 
%                                   the abf files, e.g. '20161216'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'UseOriginal': whether to use original 
%                           channel labels and units over identify_channels()
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ExpMode': experiment mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'EEG'   - EEG data; x axis in seconds; y-axis in uV
%                       'patch' - patch data; x axis in ms; y-axis in mV
%                   must be consistent with plot_traces_abf.m
%                   default == 'patch'
%                   - 'Individually': whether sweeps are plotted individually
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': the name of the directory that 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == a subdirectory named by {fileName}_traces in pwd
%                   - 'TimeUnits': units for time
%                   must be a string scalar or a character vector
%                   default == 's' for 2-data data and 'ms' for 3-data data
%                   - 'TimeStart': the start of the time interval of interest 
%                                   (in units set by TimeUnits)
%                   must be a numeric nonnegative scalar
%                   default == 0
%                   - 'TimeEnd': the end of the time interval of interest 
%                                   (in units set by TimeUnits)
%                   must be a numeric nonnegative scalar
%                   default == tVec(end)
%                   - 'ChannelTypes': the channel types
%                   must be a cellstr with nChannels elements
%                       each being one of the following:
%                           'Voltage'
%                           'Current'
%                           'Conductance'
%                           'Undefined'
%                   default == detected with identify_channels()
%                   - 'ChannelUnits': the channel units
%                   must be a cellstr with nChannels elements
%                   default == detected with identify_channels()
%                   - 'ChannelLabels': the channel labels
%                   must be a cellstr with nChannels elements
%                   default == detected with identify_channels()
%
% Requires:
%       cd/compute_and_plot_evoked_LFP.m
%       cd/compute_and_plot_concatenated_trace.m
%       cd/parse_abf.m
%       cd/plot_fields.m
%       cd/plot_traces_abf.m
%       cd/plot_FI.m
%       cd/identify_eLFP.m
%       /home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%       /home/Matlab/Downloaded_Functions/dirr.m
%       /home/Matlab/Brians_Functions/identify_CI.m
%
% File history: 
% 2016-09-22 - Created
% 2017-04-11 - Added expmode as arguments
% 2017-04-11 - Now uses dirr.m to find abf files in subdirectories too
% 2017-04-17 - BT - Creates F-I plot for current injection protocols
% 2017-04-19 - BT - Changed detection method to difference of sweep averages
% 2018-01-24 - Added isdeployed
% 2018-07-24 - Now uses a try catch statement
% 2018-09-17 - Added the input parser
% 2018-09-17 - Added identify_eLFP() and compute_and_plot_evoked_LFP()
% 2018-09-18 - Now plots a combined trace
% 2018-09-21 - Moved code to compute_and_plot_concatenated_trace.m
% 2018-09-22 - Made 'ChannelTypes', 'ChannelUnits' and 'ChannelLabels' 
%                   optional arguments
% 2018-09-22 - Added 'useOriginal' as an optional argument
% 2018-09-24 - Now tried abfload if abf2load fails

%% Hard-coded parameters
validExpModes = {'EEG', 'patch'};
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Undefined'};

%% Default values for optional arguments
directoryDefault = '';          % set later
useOriginalLabelsDefault = false;   % use identify_channels.m instead
                                    % of the original channel labels by default
expModeDefault = 'patch';       % assume traces are patching data by default
individuallyDefault = false;    % plot all sweeps together by default
outFolderDefault = '';          % set later
timeUnitsDefault = '';          % set later
timeStartDefault = [];          % set later
timeEndDefault = [];            % set later
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
if ~isdeployed
    addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));
                                            % for dirr.m, abf2load.m or abfload.m
    addpath(fullfile(functionsdirectory, '/Brians_Functions/'));        
                                            % for identify_CI.m
end

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'UseOriginal', useOriginalLabelsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ExpMode', expModeDefault, ...
    @(x) any(validatestring(x, validExpModes)));
addParameter(iP, 'Individually', individuallyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'TimeEnd', timeEndDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
useOriginal = iP.Results.UseOriginal;
expMode = validatestring(iP.Results.ExpMode, validExpModes);
individually = iP.Results.Individually;
outFolder = iP.Results.OutFolder;
timeUnits = iP.Results.TimeUnits;
timeStartUser = iP.Results.TimeStart;
timeEndUser = iP.Results.TimeEnd;
channelTypesUser = iP.Results.ChannelTypes;
channelUnitsUser = iP.Results.ChannelUnits;
channelLabelsUser = iP.Results.ChannelLabels;

% Validate channel types
if ~isempty(channelTypesUser)
    channelTypesUser = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypesUser, 'UniformOutput', false);
end

% Set dependent argument defaults
if isempty(directory)
    directory = pwd;
end

%% Find all .abf files in the directory
[~, ~, filenames] = dirr(directory, '.abf', 'name');
if isempty(filenames)
    fprintf('No abf files in current directory!\n');
    fprintf('Type ''help plot_all_abfs_dir'' for usage\n');
end
nFiles = numel(filenames);

%% Parse each file in the directory
abfParamsAllCell = cell(nFiles, 1);
dataAll = cell(nFiles, 1);
tVecAll = cell(nFiles, 1);
vVecsAll = cell(nFiles, 1);
iVecsAll = cell(nFiles, 1);
gVecsAll = cell(nFiles, 1);
dataReorderedAll = cell(nFiles, 1);
parfor iFile = 1:nFiles
%for iFile = 1:nFiles
    % Parse the abf file
    [abfParamsAllCell{iFile}, dataAll{iFile}, ...
        tVecAll{iFile}, vVecsAll{iFile}, ...
        iVecsAll{iFile}, gVecsAll{iFile}, dataReorderedAll{iFile}] = ...
        parse_abf(filenames{iFile}, 'Verbose', false, ...
                    'UseOriginal', useOriginal, ...
                    'ExpMode', expMode, ...
                    'TimeUnits', timeUnits, ...
                    'ChannelTypes', channelTypesUser, ...
                    'ChannelUnits', channelUnitsUser, ...
                    'ChannelLabels', channelLabelsUser);
end

% Convert to a struct array
abfParamsAllStruct = [abfParamsAllCell{:}];

%% Plot graphs appropriate for the identified protocol
featuresLfpAll = cell(nFiles, 1);
parfor iFile = 1:nFiles
%for iFile = 1:nFiles
%    try 
        % Extract from cell arrays
        abfParams = abfParamsAllCell{iFile};
        data = dataAll{iFile};
        iVecs = iVecsAll{iFile};
        fileName = filenames{iFile};

        % Extract some parameters
        siUs = abfParams.siUs;
        expMode = abfParams.expMode;
        timeUnits = abfParams.timeUnits;
        channelTypes = abfParams.channelTypes;
        channelUnits = abfParams.channelUnits;
        channelLabels = abfParams.channelLabels;

        % Identify whether this is a current injection protocol
        %   If so, detect spikes for each sweep and make an F-I plot
        isCI = identify_CI(iVecs, siUs);
        if isCI
            plot_FI(fileName, data, siUs);
        end

        % Identify whether this is an evoked LFP protocol
        %   If so, compute the averaged evoked LFP and plot it
        isEvokedLfp = identify_eLFP(iVecs);
        if isEvokedLfp
            % If it is an evoked LFP protocol, compute and plot it
            [tVecLfp, ~, ~, featuresLfp] = ...
                compute_and_plot_evoked_LFP(fileName, ...
                                            'OutFolder', outFolder, ...
                                            'PlotFlag', true, ...
                                            'SaveFlag', true, ...
                                            'ChannelTypes', channelTypes, ...
                                            'ChannelUnits', channelUnits, ...
                                            'ChannelLabels', channelLabels);
        else
            % Make outputs empty
            tVecLfp = [];
            featuresLfp = [];
        end

        % Set the time endpoints for individual traces
        if isEvokedLfp
            timeStart = min(tVecLfp);
            timeEnd = max(tVecLfp);
        else
            timeStart = timeStartUser;
            timeEnd = timeEndUser;
        end

        % Plot individual traces
        if ~isCI
            plot_traces_abf(fileName, 'Verbose', false, ...
                'ExpMode', expMode, 'Individually', individually, ...
                'OutFolder', outFolder, 'TimeUnits', timeUnits, ...
                'TimeStart', timeStart, 'TimeEnd', timeEnd, ...
                'ChannelTypes', channelTypes, ...
                'ChannelUnits', channelUnits, ...
                'ChannelLabels', channelLabels);            
        end

        % Save in cell arrays
        featuresLfpAll{iFile} = featuresLfp;
%   catch ME
%       fprintf('Traces for %s cannot be plotted!\n', filenames{iFile});
%       fprintf([ME.identifier, ': ', ME.message]);
%   end
end

%% Compute and plot concatenated traces for each channel
%       if the number of channels and channel types are all the same
compute_and_plot_concatenated_trace(abfParamsAllStruct, dataReorderedAll, ...
                                    'SourceDirectory', directory, ...
                                    'OutFolder', outFolder);

%% If any LFPs were computed, plot a time series for the features
% Remove empty entries
isEmpty = cellfun(@isempty, featuresLfpAll);
lfpFeaturesCell = featuresLfpAll(~isEmpty);
lfpFileNames = filenames(~isEmpty);

% Convert to a structure array
lfpFeaturesStruct = [lfpFeaturesCell{:}];

% Plot each field of the structure as its own time series
if ~isempty(lfpFeaturesStruct)
%    plot_fields(lfpFeaturesStruct, 'XTickLabel', lfpFileNames);
end

%% Copy similar figure types to its own directory
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

files = dir(directory);
if strfind(files.name, '.abf')

files = dirr(directory, '.abf');
for file = files'

    if nargin >= 3
        plot_traces_abf(fullfile(directory, filenames{iFile}), expmode, recmode);
    else
        plot_traces_abf(fullfile(directory, filenames{iFile}), expmode);
    end

    % Load abf file
    abffilename_full = construct_abffilename(filenames{iFile});    % creates full path to abf file robustly
    if exist('abf2load', 'file') == 2
        [data, siUs] = abf2load(abffilename_full);
    elseif exist('abfload', 'file') == 2
        [data, siUs] = abfload(abffilename_full);
    end

        injection_data = current_data(12000:20000,:,:);
        if std(injection_data,0,1) < 4 & abs(max(max(injection_data)) - min(min(injection_data))) > 100
            plot_FI_bt(filenames{iFile}, data, siUs);
        end
    if length(current_data) > 20000            % tests sweeps if values within injection time range are consistent
        
        injection_data = squeeze(current_data(12000:20000, 1, :));    % current values within typical injection time range
        avgs_byswp = mean(injection_data, 1);                % average of sweeps
        reduction = abs(diff(diff(avgs_byswp)));            % reduces sweeps into differences between successive sweep averages        %%% TODO: this is not differences but rather differences of differences
        [~, max_swp] = max(avgs_byswp);                    % highest sweep by average
        max_swp_peaks_avg = mean(findpeaks(injection_data(:, max_swp)));    % average peak value of greatest sweep
        if reduction < maxSwpSpacing & max_swp_peaks_avg > 100 & size(data, 3) > 1    % sweep avgs should be separated by constant
            plot_FI(filenames{iFile}, data, siUs);
        end
    end

maxSwpSpacing = 2;

function plot_all_abfs_dir (directory, expmode)

[vVecCombined] = combine_sweeps('Vectors', vVecsAll);
[iVecCombined] = combine_sweeps('Vectors', iVecsAll);
[gVecCombined] = combine_sweeps('Vectors', gVecsAll);

saveas(h, figName, 'png');

if isempty(outFolder)
    outFolder = directory;
end


%% Combine and plot all traces if the sampling intervals, channel types
%       and number of sweeps are all the same
% TODO: Put this in a function
if numel(unique(siSecondsAll)) == 1 && ...
    isequal(channelTypesAll{:}) == 1 && ...
    numel(unique(nSweepsAll)) == 1

end

vVecsAll = cell(nFiles, 1);
iVecsAll = cell(nFiles, 1);
gVecsAll = cell(nFiles, 1);
nSweepsAll = zeros(nFiles, 1);

% Combine sweeps to combine later
if nSweeps > 1
    % Decide on new nSamples and new nSweeps
    nSamplesNew = nSamples * nSweeps;
    nSweepsNew = 1;

    % Initiate new dimensions as old dimensions
    %   Note: may be 3D or 2D
    dimVVecsNew = size(vVecs);
    dimIVecsNew = size(iVecs);
    dimVVecsNew(1) = nSamplesNew;
    dimVVecsNew(2) = nSweepsNew;
    dimIVecsNew(1) = nSamplesNew;
    dimIVecsNew(2) = nSweepsNew;
    vVecsNew = squeeze(reshape(vVecs, dimVVecsNew));
    iVecsNew = squeeze(reshape(iVecs, dimIVecsNew));
    gVecsNew = [];
else
    nSamplesNew = nSamples;
    nSweepsNew = nSweeps;
    vVecsNew = vVecs;
    iVecsNew = iVecs;
    gVecsNew = gVecs;
end

vVecsAll{iFile} = vVecsNew;
iVecsAll{iFile} = iVecsNew;
gVecsAll{iFile} = gVecsNew;
nSweepsAll(iFile) = nSweepsNew;

nChannelsAll = zeros(nFiles, 1);
parfor iFile = 1:nFiles
    % Extract from arrays
    abfParams = abfParamsAll(iFile);

    % Extract from structure
    nChannelsAll{iFile} = abfParams.nChannels;
end

% Extract the common siUs
siSeconds = siSecondsAll(1);

% Count the total number of samples
nSamplesCombined = sum(nSamplesAll);

% Generate a combined time vector
tVecCombined = (1:nSamplesCombined)' * siSeconds;

% Count the maximum time
timeLength = nSamplesCombined * siSeconds;

if isempty(outputLabel)
    if ~isempty(expLabel)
        outputLabel = expLabel;
    else
        outputLabel = directoryName;
    end
end

% Combine the vectors
vVecCombined = vertcat(vVecsAll{:});
iVecCombined = vertcat(iVecsAll{:});
gVecCombined = vertcat(gVecsAll{:});

% Place vectors in cell arrays individually 
allVecs = {vVecCombined, iVecCombined, gVecCombined};
allLabels = {'Voltage (mV)', 'Current (pA)', 'Conductance (nS)'};

% Determine which vectors are not empty
nonempty = cellfun(@(x) ~isempty(x), allVecs);
nToPlot = sum(nonempty);
indToPlot = find(nonempty);

% Plot the combined trace
h = figure('Visible', 'off', 'PaperPosition', [0, 0, 80, 5]);
clf(h)
for iPlot = 1:nToPlot
    idxToPlot = indToPlot(iPlot);
    ax(iPlot) = subplot(nToPlot, 1, iPlot);
    plot(tVecCombined, allVecs{idxToPlot});
    ylabel(allLabels{idxToPlot});
    if iPlot == nToPlot
        xlabel('Time (seconds)');
    end
    xlim([min(tVecCombined), max(tVecCombined)]);
end

% Plot all traces within the LFP time window
plot_traces_abf(fileName, 'Verbose', false, ...
    'ExpMode', expMode, 'Individually', individually, ...
    'OutFolder', outFolder, 'TimeUnits', timeUnits, ...
    'TimeStart', min(tVecLfp), 'TimeEnd', max(tVecLfp), ...
    'ChannelTypes', channelTypes, ...
    'ChannelUnits', channelUnits, ...
    'ChannelLabels', channelLabels);

% If both the above are true, give preference to 
%   an evoked LFP protocol
if isCI && isEvokedLfp
    isCI = false;
end

% Save in graph outputs
graphOutputs.tVecLfp = tVecLfp;
graphOutputs.vVecLfp = vVecLfp;
graphOutputs.iVecStim = iVecStim;
graphOutputs.features = features;
graphOutputs.timeStart = timeStart;
graphOutputs.timeEnd = timeEnd;

% Save in graph outputs cell array
graphOutputsAll{iFile} = graphOutputs;

%}
