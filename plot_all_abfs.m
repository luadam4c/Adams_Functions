function [abfParamsTable, abfDataTable, abfParamsStruct, ...
            abfDataStruct, abfParamsCell, abfDataCell, ...
            featuresFileTableGabab, featuresSweepTableGabab, ...
            featuresFileTableLfp, featuresSweepTableLfp] = ...
                plot_all_abfs (varargin)
%% Plots all abf files in a directory
% Usage: [abfParamsTable, abfDataTable, abfParamsStruct, ...
%           abfDataStruct, abfParamsCell, abfDataCell] = ...
%               plot_all_abfs (varargin)
%
% Arguments:
%       varargin    - 'Directory': the name of the directory containing 
%                                   the abf files, e.g. '20161216'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileNames': names of .abf files to detect
%                   must be a cell array of character arrays or strings
%                   default == detect from pwd
%                   - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
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
%                   - 'PlotMode': plotting mode for multiple traces
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'overlapped'    - overlapped in a single plot
%                       'parallel'      - in parallel in subplots
%                   must be consistent with plot_traces_abf.m
%                   default == set in plot_traces_abf.m
%                   - 'Individually': whether sweeps are plotted individually
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'IdentifyProtocols': whether to identify protocols
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': the name of the directory that 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == same as directory
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
%                           'Other'
%                   default == detected with identify_channels()
%                   - 'ChannelUnits': the channel units
%                   must be a cellstr with nChannels elements
%                   default == detected with identify_channels()
%                   - 'ChannelLabels': the channel labels
%                   must be a cellstr with nChannels elements
%                   default == detected with identify_channels()
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for 
%                       the plot_traces() function
%
% Requires:
%       cd/all_files.m
%       cd/force_column_cell.m
%       cd/compute_and_plot_concatenated_trace.m
%       cd/parse_all_abfs.m
%       cd/plot_repetitive_protocols.m
%       cd/plot_traces_abf.m
%       cd/parse_current_injection_protocol.m
%       cd/isfigtype.m
%       /home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%
% Used by:
%       /home/zhongxiao/SCIDmiceLTP/Code/analyze_SCIDmiceLTP.m
%       /media/ashleyX/Recordings/analyze_recordings.m
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
% 2018-09-25 - Added figTypes as an argument
% 2018-09-27 - Pulled all the parsing to parse_all_abfs.m
% 2018-09-30 - Now defaults outFolder to directory
%                   but plots individual traces to subdirectories
% 2018-09-30 - Changed LFP tuning curves outFolder to fullfile(outFolder, 'LFP')
% 2018-10-03 - Updated usage of parse_all_abfs.m
%               changed outputs to allParsedParamsTable, allParsedDataTable, 
%                   abfParamsCell
% 2018-12-15 - Added 'Verbose' as a parameter
% 2018-12-15 - Fixed the case of having only one .abf file
% 2019-01-04 - Now omits the .abf files that cannot be read
% 2019-04-02 - Now does not identify protocols by default
% TODO: Restructure code so that each type of plot is its own subfunction

%% Hard-coded parameters
validExpModes = {'EEG', 'patch', ''};
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};
validPlotModes = {'overlapped', 'parallel', ''};

%% Default values for optional arguments
directoryDefault = pwd;         % look for .abf files in 
                                %   the present working directory by default
fileNamesDefault = {};          % detect from pwd by default
verboseDefault = false;         % don't print to standard output by default
useOriginalDefault = false;     % use identify_channels.m instead
                                % of the original channel labels by default
expModeDefault = 'patch';       % assume traces are patching data by default
plotModeDefault = '';           % set later
individuallyDefault = false;    % plot all sweeps together by default
identifyProtocolsDefault = false;% identify protocols by default
outFolderDefault = '';          % set later
timeUnitsDefault = '';          % set later
timeStartDefault = [];          % set later
timeEndDefault = [];            % set later
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later
figTypesDefault = 'png';
% figTypesDefault = {'png', 'fig'};

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
                                            % for abf2load.m or abfload.m
end

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'UseOriginal', useOriginalDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ExpMode', expModeDefault, ...
    @(x) isempty(x) || any(validatestring(x, validExpModes)));
addParameter(iP, 'PlotMode', plotModeDefault, ...
    @(x) isempty(x) || any(validatestring(x, validPlotModes)));
addParameter(iP, 'Individually', individuallyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'IdentifyProtocols', identifyProtocolsDefault, ...
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
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
verbose = iP.Results.Verbose;
useOriginal = iP.Results.UseOriginal;
expMode = validatestring(iP.Results.ExpMode, validExpModes);
plotMode = validatestring(iP.Results.PlotMode, validPlotModes);
individually = iP.Results.Individually;
identifyProtocols = iP.Results.IdentifyProtocols;
outFolder = iP.Results.OutFolder;
timeUnits = iP.Results.TimeUnits;
timeStartUser = iP.Results.TimeStart;
timeEndUser = iP.Results.TimeEnd;
channelTypesUser = iP.Results.ChannelTypes;
channelUnitsUser = iP.Results.ChannelUnits;
channelLabelsUser = iP.Results.ChannelLabels;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_traces() function
otherArguments = struct2arglist(iP.Unmatched);

% Validate channel types
if ~isempty(channelTypesUser)
    channelTypesUser = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypesUser, 'UniformOutput', false);
end

% Set default output folder
if isempty(outFolder)
    outFolder = directory;
end

%% Check if needed output directories exist
check_dir(outFolder);

%% Get file names
% Decide on the files to use
if isempty(fileNames)
    % Find all .abf files in the directory
    [~, fileNames] = all_files('Directory', directory, ...
                                'Extension', '.abf', ...
                                'Verbose', verbose);

    % Find all .abf files in the directory
    if isempty(fileNames)
        fprintf('Type ''help %s'' for usage\n', mfilename);
        abfParamsTable = table;
        abfDataTable = table;
        abfParamsStruct = struct;
        abfDataStruct = struct;
        abfParamsCell = {};
        abfDataCell = {};
        return
    end
end

% Force as a cell array
fileNames = force_column_cell(fileNames);

%% Parse and identify protocols from each file in the directory
% Parse all .abf files in this directory
[abfParamsTable, abfDataTable, abfParamsStruct, abfDataStruct, ...
    abfParamsCell, abfDataCell] = ...
    parse_all_abfs('FileNames', fileNames, ...
                    'Directory', directory, 'OutFolder', outFolder, ...
                    'Verbose', false, 'UseOriginal', useOriginal, ...
                    'ExpMode', expMode, 'TimeUnits', timeUnits, ...
                    'ChannelTypes', channelTypesUser, ...
                    'ChannelUnits', channelUnitsUser, ...
                    'ChannelLabels', channelLabelsUser, ...
                    'IdentifyProtocols', identifyProtocols);

% Extract each column from the data table as cell arrays
fileNames = abfParamsTable.abfFullFileName;
dataAll = abfDataTable.data;
tVecAll = abfDataTable.tVec;
vVecsAll = abfDataTable.vVecs;
iVecsAll = abfDataTable.iVecs;
gVecsAll = abfDataTable.gVecs;
dataReorderedAll = abfDataTable.dataReordered;

% Force as a cell array again
fileNames = force_column_cell(fileNames);

% Count the number of files
nFiles = numel(fileNames);

%% Plot F-I plots
if identifyProtocols
    parfor iFile = 1:nFiles
        % Extract from cell arrays
        abfParams = abfParamsCell{iFile};

        % Extract some parameters
        isCI = abfParams.isCI;

        % Only do anything for current injection protocols
        if isCI
            % Extract from cell arrays
            data = dataAll{iFile};
            fileName = fileNames{iFile};

            % Extract some more parameters
            siUs = abfParams.siUs;

            % Detect spikes for each sweep and make an F-I plot
            parse_current_injection_protocol(fileName, data, siUs);
        end
        
        close all force hidden;
    end
end

%% Deal with special protocols
if identifyProtocols
    % Plot all graphs associated with GABA-B IPSC protocols
    [featuresFileTableGabab, featuresSweepTableGabab] = ...
        plot_repetitive_protocols('EvokedGABAB', ...
                        'FileNames', fileNames, 'OutFolder', outFolder, ...
                        'Verbose', verbose, 'UseOriginal', useOriginal, ...
                        'ExpMode', expMode, 'PlotMode', plotMode, ...
                        'Individually', individually, 'OutFolder', outFolder, ...
                        'TimeUnits', timeUnits, 'FigTypes', figTypes, ...
                        'ChannelTypes', channelTypesUser, ...
                        'ChannelUnits', channelUnitsUser, ...
                        'ChannelLabels', channelLabelsUser, ...
                        'AbfParamsCell', abfParamsCell, 'AbfDataCell', abfDataCell);

    % Plot all graphs associated with LFP protocols
    [featuresFileTableLfp, featuresSweepTableLfp] = ...
        plot_repetitive_protocols('EvokedLFP', ...
                        'FileNames', fileNames, 'OutFolder', outFolder, ...
                        'Verbose', verbose, 'UseOriginal', useOriginal, ...
                        'ExpMode', expMode, 'PlotMode', plotMode, ...
                        'Individually', individually, 'OutFolder', outFolder, ...
                        'TimeUnits', timeUnits, 'FigTypes', figTypes, ...
                        'ChannelTypes', channelTypesUser, ...
                        'ChannelUnits', channelUnitsUser, ...
                        'ChannelLabels', channelLabelsUser, ...
                        'AbfParamsCell', abfParamsCell, 'AbfDataCell', abfDataCell);
end

%% Plot individual traces
%parfor iFile = 1:nFiles
for iFile = 1:nFiles
    % Extract from cell arrays
    abfParams = abfParamsCell{iFile};
    abfData = abfDataCell{iFile};

    % Extract some parameters
    if identifyProtocols
        isCI = abfParams.isCI;
        isEvokedLfp = abfParams.isEvokedLfp;
        isEvokedGabab = abfParams.isEvokedGabab;
    end
    
    % Don't do this for current injection protocols
    %   or evoked GABA-B IPSC protocols or evoked LFP protocols
    %   because individual traces are already plotted
    if ~identifyProtocols || (~isCI && ~isEvokedLfp && ~isEvokedGabab)
        % Extract from cell arrays
        fileName = fileNames{iFile};

        % Extract some more parameters
        expMode = abfParams.expMode;
        timeUnits = abfParams.timeUnits;
        channelTypes = abfParams.channelTypes;
        channelUnits = abfParams.channelUnits;
        channelLabels = abfParams.channelLabels;

        % Set the time endpoints for individual traces
        timeStart = timeStartUser;
        timeEnd = timeEndUser;

        % Plot individual traces
        %   Note: Must make outFolder empty so that outputs
        %           be plotted in subdirectories
        plot_traces_abf(fileName, 'Verbose', false, ...
            'ParsedParams', abfParams, 'ParsedData', abfData, ...
            'ExpMode', expMode, 'PlotMode', plotMode, ...
            'Individually', individually, ...
            'OutFolder', '', 'TimeUnits', timeUnits, ...
            'TimeStart', timeStart, 'TimeEnd', timeEnd, ...
            'ChannelTypes', channelTypes, ...
            'ChannelUnits', channelUnits, ...
            'ChannelLabels', channelLabels, ...
            'FigTypes', figTypes, otherArguments{:});
    end
    
    close all force hidden;
end

%% Compute and plot concatenated traces for each channel
%       if the number of channels and channel types are all the same

compute_and_plot_concatenated_trace(abfParamsStruct, dataReorderedAll, ...
                                    'SourceDirectory', directory, ...
                                    'OutFolder', outFolder, otherArguments{:});

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
        plot_traces_abf(fullfile(directory, fileNames{iFile}), expmode, recmode);
    else
        plot_traces_abf(fullfile(directory, fileNames{iFile}), expmode);
    end

    % Load abf file
    abffilename_full = construct_abffilename(fileNames{iFile});    % creates full path to abf file robustly
    if exist('abf2load', 'file') == 2
        [data, siUs] = abf2load(abffilename_full);
    elseif exist('abfload', 'file') == 2
        [data, siUs] = abfload(abffilename_full);
    end

        injection_data = current_data(12000:20000,:,:);
        if std(injection_data,0,1) < 4 & abs(max(max(injection_data)) - min(min(injection_data))) > 100
            parse_current_injection_protocol_bt(fileNames{iFile}, data, siUs);
        end
    if length(current_data) > 20000            % tests sweeps if values within injection time range are consistent
        
        injection_data = squeeze(current_data(12000:20000, 1, :));    % current values within typical injection time range
        avgs_byswp = mean(injection_data, 1);                % average of sweeps
        reduction = abs(diff(diff(avgs_byswp)));            % reduces sweeps into differences between successive sweep averages        %%% TODO: this is not differences but rather differences of differences
        [~, max_swp] = max(avgs_byswp);                    % highest sweep by average
        max_swp_peaks_avg = mean(findpeaks(injection_data(:, max_swp)));    % average peak value of greatest sweep
        if reduction < maxSwpSpacing & max_swp_peaks_avg > 100 & size(data, 3) > 1    % sweep avgs should be separated by constant
            parse_current_injection_protocol(fileNames{iFile}, data, siUs);
        end
    end

maxSwpSpacing = 2;

function plot_all_abfs (directory, expmode)

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

%% Plot graphs appropriate for the identified protocol
parfor iFile = 1:nFiles
   try
   catch ME
       fprintf('Traces for %s cannot be plotted!\n', fileNames{iFile});
       fprintf([ME.identifier, ': ', ME.message]);
   end
end

%       /home/Matlab/Brians_Functions/identify_CI.m

% Set dependent argument defaults
if isempty(directory)
    directory = pwd;
end

%% Identify protocols
isCIAll = false(nFiles, 1);
isEvokedLfpAll = false(nFiles, 1);
parfor iFile = 1:nFiles
        % Extract from cell arrays
        abfParams = abfParamsCell{iFile};
        iVecs = iVecsAll{iFile};

        % Extract some parameters
        siUs = abfParams.siUs;

        % Identify whether this is a current injection protocol
        isCIAll(iFile) = identify_CI(iVecs, siUs);

        % Identify whether this is an evoked LFP protocol
        isEvokedLfpAll(iFile) = identify_eLFP(iVecs);
end

% Decide on the files to use
if isempty(fileNames)
    % Find all .abf files in the directory
    [~, ~, fileNames] = dirr(directory, '.abf', 'name');
    if isempty(fileNames)
        fprintf('No abf files in current directory!\n');
        fprintf('Type ''help plot_all_abfs'' for usage\n');
        return
    end
end

fileNames = cellfun(@(x, y) fullfile(x, y), ...
                    {files.folder}, {files.name}, ...
                    'UniformOutput', false);

lfpFigNames = ...
    cellfun(@(x) fullfile(outFolder, [x, '_vs_', xLabel]), ...
            lfpFieldNames, 'UniformOutput', false);


%       cd/get_column.m
dataAll = get_column(abfDataTable, 'data');
tVecAll = get_column(abfDataTable, 'tVec');
vVecsAll = get_column(abfDataTable, 'vVecs');
iVecsAll = get_column(abfDataTable, 'iVecs');
gVecsAll = get_column(abfDataTable, 'gVecs');
dataReorderedAll = get_column(abfDataTable, 'dataReordered');

%       cd/extract_fullpath.m

if isempty(fileNames)
    % Find all .abf files in the directory
    files = dir(fullfile(directory, '*.abf'));
    if isempty(files)
        fprintf('No .abf files in current directory!\n');
        fprintf('Type ''help %s'' for usage\n', mfilename);
        abfParamsStruct = struct;
        dataAll = {};
        tVecAll = {};
        vVecsAll = {};
        iVecsAll = {};
        gVecsAll = {};
        dataReorderedAll = {};
        abfParamsCell = {};
        return
    end

    % Construct the full file names
    fileNames = extract_fullpath(files);
end

if isEvokedLfp
    timeStart = min(tVecLfp);
    timeEnd = max(tVecLfp);
else
    timeStart = timeStartUser;
    timeEnd = timeEndUser;
end

plot_struct(lfpFeaturesStruct, ...
                'XTickLabels', xTickLabels, ...
                'XLabel', xLabel, ...
                'FigNames', lfpFigNames);


% Get all field names
lfpFieldNames = fieldnames(lfpFeaturesStruct);

% Create figure names
lfpFigNames = ...
    cellfun(@(x) fullfile(outFolderLfp, [x, '_vs_', xLabel]), ...
            lfpFieldNames, 'UniformOutput', false);

% Remove empty entries
isEmpty = cellfun(@isempty, featuresLfpCell);
lfpFeaturesCell = featuresLfpCell(~isEmpty);
lfpFileNames = fileNames(~isEmpty);

if ~isempty(lfpFeaturesStruct)

outFolderLfp = fullfile(outFolder, 'LFPs');

featuresTableGabab = ...
featuresTableLfp = ...

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
