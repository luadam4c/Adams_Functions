function [featuresFileTable, featuresSweepTable] = ...
                plot_protocols (protocolType, varargin)
%% Computes features for each file according to protocol type
% Usage: [featuresFileTable, featuresSweepTable] = ...
%               plot_protocols (protocolType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       featuresFileTable   - features table for each file
%                           specified as a 2-D table
%       featuresSweepTable  - features table for each sweep
%                           specified as a 2-D table
% Arguments:
%       protocolType- protocol type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'EvokedLFP'   - an evoked LFP protocol
%                       'EvokedGABAB' - an evoked GABA-B protocol
%                   must be consistent with plot_traces_abf.m
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
%                   default == 'overlapped'
%                   - 'Individually': whether sweeps are plotted individually
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SavePlotsFlag': whether to save plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveTablesFlag': whether to save tables
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'OutFolder': the name of the directory that 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == same as directory
%                   - 'TimeUnits': units for time
%                   must be a string scalar or a character vector
%                   default == 's' for 2-data data and 'ms' for 3-data data
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
%                   - 'AbfParamsCell': parsed .abf file parameters
%                   must be a cell array of structures or tables
%                   default == loaded from fileNames
%                   - 'AbfDataCell': parsed .abf file data
%                   must be a cell array of structures or tables
%                   default == loaded from fileNames
%                   - Any other parameter-value pair for the plot() function
%                   
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/compute_and_plot_all_responses.m
%       cd/compute_and_plot_average_response.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/parse_all_abfs.m
%       cd/plot_table.m
%       cd/plot_traces_abf.m
%
% Used by:
%       cd/plot_all_abfs.m

% File History:
% 2018-12-15 Created by Adam Lu
% 2018-12-17 Added savePlotsFlag and saveTablesFlag
% 2018-12-18 Moved code to plot_table.m
% TODO: Make plot flags for each type of plot optional arguments
% TODO: Decide on whether to use parfor across files
%       based on the total number of sweeps versus the total number of files
%       cf. m3ha
% 

%% TODO: Make these parameters
plotSeparateFlag = false; %true;
plotAltogetherFlag = false; %true;
plotAverageFlag = false; %true;

%% Hard-coded parameters
validProtocolTypes = {'EvokedLFP', 'EvokedGABAB'};
validExpModes = {'EEG', 'patch', ''};
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};
validPlotModes = {'overlapped', 'parallel'};
featuresFileSheetName = 'features_table_by_file.xlsx';
featuresSweepSheetName = 'features_table_by_sweep.xlsx';
xLabelFile = 'fileName';
xLabelSweep = 'sweepNumber';

% For Evoked GABA-B IPSC
lowPassFrequencyGabab = 1000;       % lowpass filter frequency in Hz
medFiltWindowGabab = [];            % median filter window in ms
smoothWindowGabab = 100;            % moving average filter window in ms
baselineLengthMsGabab = 50;         % baseline length in ms
responseLengthMsGabab = 500;        % response length in ms
minPeakDelayMsGabab = 50;           % min peak delay after pulse end in ms

% For Evoked LFP
lowPassFrequencyLfp = 1000;         % lowpass filter frequency in Hz
medFiltWindowLfp = [];              % median filter window in ms
smoothWindowLfp = [];               % moving average filter window in ms
baselineLengthMsLfp = 50;           % baseline length in ms
responseLengthMsLfp = 20;           % response length in ms
minPeakDelayMsLfp = 1;              % min peak delay after pulse end in ms
peakDirectionLfp = 'upward';

%% Default values for optional arguments
directoryDefault = pwd;         % look for .abf files in 
                                %   the present working directory by default
fileNamesDefault = {};          % detect from pwd by default
verboseDefault = false;         % don't print to standard output by default
useOriginalDefault = false;     % use identify_channels.m instead
                                % of the original channel labels by default
expModeDefault = 'patch';       % assume traces are patching data by default
plotModeDefault = 'overlapped'; % plot traces overlapped by default
individuallyDefault = false;    % plot all sweeps together by default
savePlotsFlagDefault = true;    % save all plots by default
saveTablesFlagDefault = true;   % save all tables by default
outFolderDefault = '';          % set later
timeUnitsDefault = '';          % set later
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later
figTypesDefault = {'png'};
% figTypesDefault = {'png', 'fig'};
abfParamsCellDefault = {};      % set later
abfDataCellDefault = {};        % set later

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
addRequired(iP, 'protocolType', ...
    @(x) any(validatestring(x, validProtocolTypes)));

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
    @(x) any(validatestring(x, validExpModes)));
addParameter(iP, 'PlotMode', plotModeDefault, ...
    @(x) any(validatestring(x, validPlotModes)));
addParameter(iP, 'Individually', individuallyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SavePlotsFlag', savePlotsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveTablesFlag', saveTablesFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'AbfParamsCell', abfParamsCellDefault, ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));
addParameter(iP, 'AbfDataCell', abfDataCellDefault, ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));

% Read from the Input Parser
parse(iP, protocolType, varargin{:});
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
verbose = iP.Results.Verbose;
useOriginal = iP.Results.UseOriginal;
expMode = validatestring(iP.Results.ExpMode, validExpModes);
plotMode = validatestring(iP.Results.PlotMode, validPlotModes);
individually = iP.Results.Individually;
savePlotsFlag = iP.Results.SavePlotsFlag;
saveTablesFlag = iP.Results.SaveTablesFlag;
outFolder = iP.Results.OutFolder;
timeUnits = iP.Results.TimeUnits;
channelTypesUser = iP.Results.ChannelTypes;
channelUnitsUser = iP.Results.ChannelUnits;
channelLabelsUser = iP.Results.ChannelLabels;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
abfParamsCell = iP.Results.AbfParamsCell;
abfDataCell = iP.Results.AbfDataCell;

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

% Validate Protocol type
protocolType = validatestring(protocolType, validProtocolTypes);

% Validate channel types
if ~isempty(channelTypesUser)
    channelTypesUser = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypesUser, 'UniformOutput', false);
end

% Set default output folder
if isempty(outFolder)
    outFolder = directory;
end

%% Set parameters based on protocol type
switch protocolType
    case 'EvokedLFP'
        % General
        isProtocolStr = 'isEvokedLfp';

        % For computing
        responseType = 'Voltage';
        peakDirection = peakDirectionLfp;
        lowPassFrequency = lowPassFrequencyLfp;
        medFiltWindow = medFiltWindowLfp;
        smoothWindow = smoothWindowLfp;
        baselineLengthMs = baselineLengthMsLfp;
        responseLengthMs = responseLengthMsLfp;
        minPeakDelayMs = minPeakDelayMsLfp;

        % For plotting
        outFolderProtocolName = 'LFPs';
        varToPlot = {'peakAmplitude', 'peakDelayMs', ...
                    'peakHalfWidthMs', 'peakTimeConstantMs', 'peak90DecayMs'};
        fileSuffix = '_LFP';
        responseName = 'Evoked potential';
    case 'EvokedGABAB'
        % General
        isProtocolStr = 'isEvokedGabab';

        % For computing
        responseType = 'Current';
        peakDirection = 'upward';
        lowPassFrequency = lowPassFrequencyGabab;
        medFiltWindow = medFiltWindowGabab;
        smoothWindow = smoothWindowGabab;
        baselineLengthMs = baselineLengthMsGabab;
        responseLengthMs = responseLengthMsGabab;
        minPeakDelayMs = minPeakDelayMsGabab;

        % For plotting
        outFolderProtocolName = 'GABAB-IPSCs';
        varToPlot = {'peakAmplitude', 'peakDelayMs', ...
                    'peakHalfWidthMs', 'peakTimeConstantMs', 'peak90DecayMs'};
        fileSuffix = '_GABAB';
        responseName = 'GABA_B IPSC';
    otherwise
        error('Code logic error!!\n');
end

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
        featuresFileTable = table;
        return
    end
end

% Count the number of files
nFiles = numel(fileNames);

% Get outFolder name
[~, outFolderName] = fileparts(outFolder);

%% Parse and identify protocols from each file in the directory
% Parse all .abf files if not already done 
if isempty(abfParamsCell) || isempty(abfDataCell)
    [~, ~, ~, ~, abfParamsCell, abfDataCell] = ...
        parse_all_abfs('FileNames', fileNames, 'OutFolder', outFolder, ...
                        'Verbose', false, 'UseOriginal', useOriginal, ...
                        'ExpMode', expMode, 'TimeUnits', timeUnits, ...
                        'ChannelTypes', channelTypesUser, ...
                        'ChannelUnits', channelUnitsUser, ...
                        'ChannelLabels', channelLabelsUser, ...
                        'IdentifyProtocols', true);
end

% Set the output folder
outFolderProtocol = fullfile(outFolder, outFolderProtocolName);

%% Do the job
% Compute features and plot protocol traces
featuresPerFileCell = cell(nFiles, 1);
featuresPerSweepCell = cell(nFiles, 1);
for iFile = 1:nFiles
%parfor iFile = 1:nFiles
    % Extract from cell arrays
    abfParams = abfParamsCell{iFile};

    % Extract some parameters
    isProtocol = abfParams.(isProtocolStr);

    % Only do anything for a matching protocol
    if isProtocol
        % Extract from cell arrays
        fileName = fileNames{iFile};
        abfData = abfDataCell{iFile};

        % Extract some more parameters
        timeUnits = abfParams.timeUnits;
        channelTypes = abfParams.channelTypes;
        channelUnits = abfParams.channelUnits;
        channelLabels = abfParams.channelLabels;

        % Compute the average response, parse the features and plot
        [tVecAvg, respAvg, stimAvg, featuresAvg, h1] = ...
            compute_and_plot_average_response(fileName, responseType, ...
                'LowPassFrequency', lowPassFrequency, ...
                'MedFiltWindow', medFiltWindow, ...
                'SmoothWindow', smoothWindow, ...
                'BaselineLengthMs', baselineLengthMs, ...
                'ResponseLengthMs', responseLengthMs, ...
                'MinPeakDelayMs', minPeakDelayMs, ...
                'PeakDirection', peakDirection, ...
                'OutFolder', outFolderProtocol, ...
                'FileSuffix', fileSuffix, 'ResponseName', responseName, ...
                'PlotFlag', plotAverageFlag, 'SavePlotsFlag', savePlotsFlag, ...
                'SaveTablesFlag', saveTablesFlag, 'FigTypes', figTypes, ...
                'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
                'ChannelLabels', channelLabels, ...
                'ParsedParams', abfParams, 'ParsedData', abfData, ...
                otherArguments);

        % Parse features from and plot all responses
        [tVecAll, respAll, stimAll, featuresAll, h2] = ...
            compute_and_plot_all_responses(fileName, responseType, ...
                'LowPassFrequency', lowPassFrequency, ...
                'MedFiltWindow', medFiltWindow, ...
                'SmoothWindow', smoothWindow, ...
                'BaselineLengthMs', baselineLengthMs, ...
                'ResponseLengthMs', responseLengthMs, ...
                'MinPeakDelayMs', minPeakDelayMs, ...
                'PeakDirection', peakDirection, ...
                'FileSuffix', fileSuffix, 'ResponseName', responseName, ...
                'PlotFlag', plotSeparateFlag, 'SavePlotsFlag', savePlotsFlag, ...
                'SaveTablesFlag', saveTablesFlag, 'FigTypes', figTypes, ...
                'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
                'ChannelLabels', channelLabels, ...
                'ParsedParams', abfParams, 'ParsedData', abfData, ...
                otherArguments);

        % Plot individual protocol traces, stimulus and response separate
        %   Note: Must make outFolder empty so that outputs
        %           be plotted in subdirectories
        if plotAltogetherFlag
            % Set the time endpoints for individual protocol traces
            [timeStart, timeEnd] = ...
                argfun(@(x) extract_elements(tVecAll, x), 'first', 'last');

            % Plot individual traces
            plot_traces_abf(fileName, ...
                'ParsedParams', abfParams, 'ParsedData', abfData, ...
                'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
                'ChannelLabels', channelLabels, ...
                'Verbose', verbose, 'ExpMode', expMode, ...
                'PlotMode', plotMode, 'Individually', individually, ...
                'OutFolder', '', 'TimeUnits', timeUnits, ...
                'TimeStart', timeStart, 'TimeEnd', timeEnd, ...
                'FigTypes', figTypes);
        end
    else
        % Make outputs an empty table
        featuresAvg = [];
        featuresAll = [];
    end


    % Save feature table in cell array
    featuresPerFileCell{iFile} = featuresAvg;
    featuresPerSweepCell{iFile} = featuresAll;

    % Close figures before exiting parfor loop
    close all force hidden
end

% Remove empty entries
isEmpty = cellfun(@isempty, featuresPerFileCell);
featuresFileProtocol = featuresPerFileCell(~isEmpty);
featuresSweepProtocol = featuresPerSweepCell(~isEmpty);

% Join the tables together
featuresFileTable = vertcat(featuresFileProtocol{:});
featuresSweepTable = vertcat(featuresSweepProtocol{:});

% Create file paths for the tables
featuresFilePath = fullfile(outFolder, [outFolderName, '_', ...
                        outFolderProtocolName, '_', featuresFileSheetName]);
featuresSweepPath = fullfile(outFolder, [outFolderName, '_', ...
                        outFolderProtocolName, '_', featuresSweepSheetName]);

% If any features were computed, 
%   plot each column of the feature table as its own time series
if ~isempty(featuresFileTable) && istable(featuresFileTable)
    % Save the table
    writetable(featuresFileTable, featuresFilePath, 'WriteRowNames', true);

    % Plot each variable of the table
    plot_table(featuresFileTable, 'VariableNames', varToPlot, ...
               'OutFolder', outFolderProtocol, 'XLabel', xLabelFile);
end
if ~isempty(featuresSweepTable) && istable(featuresSweepTable)
    % Save the table
    writetable(featuresSweepTable, featuresSweepPath, 'WriteRowNames', false);

    % Plot each variable of the table
    plot_table(featuresSweepTable, 'VariableNames', varToPlot, ...
               'OutFolder', outFolderProtocol, 'XLabel', xLabelSweep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute the averaged evoked LFP and plot it
[tVecAvg, ~, ~, featuresAvg] = ...
    compute_and_plot_evoked_LFP(fileName, ...
        'ParsedParams', abfParams, 'ParsedData', abfData, ...                
        'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
        'ChannelLabels', channelLabels, 'OutFolder', outFolderProtocol, ...
        'PlotFlag', true, 'SavePlotsFlag', true);

% Use file names as row names
featuresFileTable.Properties.RowNames = fileNamesProtocol;

% Set the time endpoints for protocol traces, averaged or not
timeStart = min(tVecAll);
timeEnd = max(tVecAll);

fileNamesProtocol = fileNames(~isEmpty);

nSweepsPreviousFiles = 0;   % counts the number of sweeps in previous files
'MinRowNumber', nSweepsPreviousFiles + 1, ...
% Update number of sweeps in previous files
nSweepsPreviousFiles = nSweepsPreviousFiles + nVectors;

timeStart = cellfun(@min, tVecAll);
timeEnd = cellfun(@max, tVecAll);

nVectors = count_vectors(tVecAll);

params.FileBase = [fileBase, '_Swp', num2str(iVec)];

timeStart = extract_elements(tVecAll, 'first');
timeEnd = extract_elements(tVecAll, 'last');

% Save in a single params structure
params = table2struct(featuresAll(iVec, :));
params.OutFolder = outFolder;
params.SaveFlag = savePlotsFlag;
params.FigTypes = 'png';
params.FileBase = newFileBases{iVec};
params.FileSuffix = fileSuffix;
params.ResponseName = responseName;

plot_table(featuresFileTable, varToPlot, outFolderProtocol, xLabelFile);
plot_table(featuresSweepTable, varToPlot, outFolderProtocol, xLabelSweep);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
