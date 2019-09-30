function handles = plot_traces_spike2_mat (spike2Path, varargin)
%% Plots traces from a Spike2-exported .mat file
% Usage: handles = plot_traces_spike2_mat (spike2Path, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [~, matPaths] = all_files('Ext', 'mat');
%       for i = 1:numel(matPaths); plot_traces_spike2_mat(matPaths{i}, 'ParseLaser', true); end
%
% Outputs:
%       handles     - graphics object handles, containing:
%                       fig
%                       ax
%                   specified as a structure
%
% Arguments:
%       spike2Path  - path to Spike2-exported .mat file
%                   must be a string scalar or a character vector
%       varargin    - 'ParseGas': whether to parse pleth pulses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ParseLaser': whether to parse laser pulses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TimeUnits': output time units
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'min'   - minutes
%                       's'     - seconds
%                       'ms'    - milliseconds
%                       'us'    - microseconds
%                   default == 's'
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == {'png', 'epsc'}
%                   - Any other parameter-value pair for plot_traces() TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_subplots.m
%       cd/create_time_vectors.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/find_in_strings.m
%       cd/find_matching_files.m
%       cd/find_window_endpoints.m
%       cd/force_matrix.m
%       cd/force_string_end.m
%       cd/parse_spike2_mat.m
%       cd/plot_vertical_line.m
%       cd/plot_window_boundaries.m
%       cd/set_figure_properties.m
%       cd/save_all_figtypes.m
%       cd/update_figure_for_corel.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-30 Moved from plethRO1_analyze.m
% 

%% Hard-coded constants
S_PER_MIN = 60;

%% Hard-coded parameters
validTimeUnits = {'min', 's', 'ms', 'us'};

% TODO: Make optional arguments
channelNamesUser = {'Sound'; 'Pleth 2'; 'WIC#2'; 'WIC#1'};
plotChannelNames = {'Laser'; 'Pleth'; 'EEG'; 'EMG'};
plotSpectrogram = true;
plotZooms = false;
alignToFirstStim = true; %false;
plotForCorel = false; %true;
removeTicks = false;
outFolder = pwd;
figBase = '';
figTypes = {'png', 'epsc2'};
stimTableSuffix = '_pulses';
spectYLimits = [1, 50];
relativeWindow = [-10, 10];
% relativeWindowNoSwd = ;
% relativeWindowSwd = ;
% relativeWindow = [-10, 20];
% relativeWindowNoSwd = [-6.3, -6];
% relativeWindowSwd = [6, 6.3];
suffixNoSwd = '_no_swd';
suffixSwd = '_swd';

%% Default values for optional arguments
parseGasDefault = false;
parseLaserDefault = false;
timeUnitsDefault = '';          % set later
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'spike2Path', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ParseGas', parseGasDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ParseLaser', parseLaserDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, spike2Path, varargin{:});
parseGas = iP.Results.ParseGas;
parseLaser = iP.Results.ParseLaser;
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_traces() function
% otherArguments = iP.Unmatched;

%% Preparation
% Force the Spike2 path to end with '.mat'
spike2Path = force_string_end(spike2Path, '.mat');

% Decide on a figure base
if isempty(figBase)
    figBase = extract_fileparts(spike2Path, 'base');
end

% Decide on the time units
if parseGas
    timeUnits = 'min';
elseif parseLaser
    timeUnits = 's';
else
    timeUnits = 's';
end

% Compute the time units conversion factor
switch timeUnits
    case 's'
        unitsPerSecond = 1;
    case 'min'
        unitsPerSecond = 1 / S_PER_MIN;
    otherwise
        error('timeUnits unrecognized!');
end

% Decide on the figure expansion
if plotForCorel
    figExpansion = [1, 0.4];
else
    figExpansion = [1, 1];
end

% Decide on a figure base
figBaseNoSwd = strcat(figBase, suffixNoSwd);
figBaseSwd = strcat(figBase, suffixSwd);

% Create figure paths
[figPathBase, figPathBaseNoSwd, figPathBaseSwd] = ...
    argfun(@(x) fullfile(outFolder, x), figBase, figBaseNoSwd, figBaseSwd);

%% Load trace data
% Load data and parse gas trace if necessary
data = parse_spike2_mat(spike2Path, 'ChannelNames', channelNamesUser, ...
                        'ParseGas', parseGas, 'ParseLaser', parseLaser);

% Extract the channel values and force as a matrix
channelValues = force_matrix(data.channelValues);

% Find the common sampling interval in seconds
siSeconds = nanmean(data.siSeconds);

% Find the common number of samples
nSamples = nanmean(data.nSamples);

% Create a time vector in minutes
timeVec = create_time_vectors(nSamples, 'TimeUnits', timeUnits, ...
                                'SamplingIntervalSec', siSeconds);

%% Restrict or align trace data if requested
if alignToFirstStim
    % Find corresponding stim table
    [~, stimTablePath] = ...
        find_matching_files(spike2Path, 'Suffix', stimTableSuffix, ...
                            'Extension', 'csv');

    % Extract the first stim start time from the stim table
    stimTable = readtable(stimTablePath);
    startTimesSec = stimTable.startTime;
    startTimes = startTimesSec * unitsPerSecond;
    centerTime = startTimes(1);

    % Create time window
    timeWindow = centerTime + relativeWindow;

    % Find the time window endpoints
    endPoints = find_window_endpoints(timeWindow, timeVec);

    % Restrict time vector and shift it to be 
    %   relative to the first stim start time
    timeVecToPlot = extract_subvectors(timeVec, 'EndPoints', endPoints);
    timeVecToPlot = timeVecToPlot - centerTime;

    % Restrict data to the time window endpoints
    channelValuesToPlot = ...
        extract_subvectors(channelValues, 'EndPoints', endPoints);

else
    timeVecToPlot = timeVec;
    channelValuesToPlot = channelValues;
end

%% Compute the spectrogram from the EEG trace
% Look for an EEG trace
eegChannelNum = find_in_strings('EEG', plotChannelNames, 'MaxNum', 1);

% Only compute spectrogram if EEG trace is found
if plotSpectrogram && ~isempty(eegChannelNum)
    % Update channel names to plot
    plotChannelNames = [plotChannelNames; {'Spect'}];

    % Get the EEG trace
    eegToPlot = channelValuesToPlot(:, eegChannelNum);

    % Compute the spectrogram
    [spectData, freqHz, timeInstantsRelSeconds] = ...
        compute_spectrogram(eegToPlot, siSeconds);

    % Convert times to minutes
    timeInstantsRel = timeInstantsRelSeconds * unitsPerSecond;

    % Align times if requested
    if alignToFirstStim
        timeInstantsMin = relativeWindow(1) + timeInstantsRel;
    else
        timeInstantsMin = timeInstantsRel;
    end
else
    plotSpectrogram = false;
end

%% Plot traces
% Count the number of subplots
nSubPlots = numel(plotChannelNames);

% Create figure
[fig, ax] = create_subplots(nSubPlots, 1, 'AlwaysNew', true, ...
                            'FigExpansion', figExpansion);

% TODO: Let plot_traces accept 'AxesHandles' and use it
% Plot traces
for iPlot = 1:nSubPlots
    % Plot the appropriate trace or map
    if plotSpectrogram && iPlot == nSubPlots
        % TODO: plot_spectrogram(spectData, timeInstantsMin, freqHz, 'AxesHandle', ax);
        spectColorMapFile = '/media/adamX/Settings_Matlab/spectrogram_colormap.mat';

        % Plot the log spectrum
        imagesc(ax(iPlot), timeInstantsMin, freqHz, abs(spectData));

        % Set a colormap
        colorMapSpectFile = matfile(spectColorMapFile);
        colorMapSpect = colorMapSpectFile.colorMap;
        colormap(ax(iPlot), colorMapSpect);

        % Flip the Y Axis so lower frequencies are at the bottom
        set(ax(iPlot), 'YDir', 'normal');

        % Restrict to certain y axis limits
        set(ax(iPlot), 'YLim', spectYLimits);
    else
        plot(ax(iPlot), timeVecToPlot, channelValuesToPlot(:, iPlot), ...
                'Color', 'k');
        if iPlot == 1
%            set(ax(iPlot), 'YLim', [33, 40]);
        elseif iPlot == 2
%            set(ax(iPlot), 'YLim', [-1, 1]);
        elseif iPlot == 3
            % Show only -0.5-0.5 mV
%            set(ax(iPlot), 'YLim', [-0.5, 0.5]);
        elseif iPlot == 4
%            set(ax(iPlot), 'YLim', [-0.5, 0.5]);
        end
    end

    % Plot a vertical line
    if alignToFirstStim
        plot_vertical_line(0, 'AxesHandle', ax(iPlot), 'Color', 'k', ...
                            'LineWidth', 1);
    end
end

% Link the x axes
linkaxes(ax, 'x');

% Update figure for Corel Draw
if plotForCorel
    fig = update_figure_for_corel(fig, 'RemoveTicks', removeTicks);
end

% Save the figure in different zooms
if plotZooms
    % Shrink figure width
    set_figure_properties('FigHandle', fig, 'FigExpansion', [0.5, 1], ...
                            'ExpandFromDefault', false);

    % Zoom to first region and save
    xlim(relativeWindowNoSwd);
    save_all_figtypes(fig, figPathBaseNoSwd, figTypes);

    % Zoom to second region and save
    xlim(relativeWindowSwd);
    save_all_figtypes(fig, figPathBaseSwd, figTypes);

    % Expand figure width
    set_figure_properties('FigHandle', fig, 'FigExpansion', [2, 1], ...
                            'ExpandFromDefault', false);

    % Plot shades
    for iPlot = 1:numel(plotChannelNames)
        % Plot vertical shades
        subplot(ax(iPlot))
        plot_window_boundaries([relativeWindowNoSwd, relativeWindowSwd], ...
                                'BoundaryType', 'verticalShade');
    end

    % Zoom back out
    xlim(relativeWindow);
end

% Save figure
save_all_figtypes(fig, figPathBase, figTypes);    

%% Output results
handles.fig = fig;
handles.ax = ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spectData, freqHz, timeInstantsSeconds] = ...
                compute_spectrogram(eegToPlot, siSeconds)
%% Computes a spectrogram
% TODO: Pull out as its own function

%% Hard-coded parameters
% TODO: Make optional arguments
binWidthSeconds = 1;             % 1 second windows
overlapSeconds = [];

%% Preparation
% Compute the bin width in samples
binWidthSamples = round(binWidthSeconds / siSeconds);

% Set a default overlap in samples
if isempty(overlapSeconds)
    overlapSamples = round(binWidthSamples / 2);
else
    overlapSamples = round(overlapSeconds / siSeconds);
end

% Compute the sampling frequency in Hz
samplingFreqHz = 1 / siSeconds;

%% Do the job
% Compute the spectrogram
[spectData, freqHz, timeInstantsSeconds] = ...
    spectrogram(eegToPlot, binWidthSamples, ...
                overlapSamples, [], samplingFreqHz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

spike2PathBase = extract_fileparts(spike2Path, 'pathbase');
stimTablePath = [spike2PathBase, stimTableSuffix, '.csv'];

% Extract the channel names and values
channelNames = data.channelNames;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%