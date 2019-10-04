function handles = plot_traces_spike2_mat (spike2Path, varargin)
%% Plots traces from a Spike2-exported .mat file
% Usage: handles = plot_traces_spike2_mat (spike2Path, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [~, matPaths] = all_files('Ext', 'mat');
%       for i = 1:numel(matPaths); plot_traces_spike2_mat(matPaths{i}, 'AlignToStim', true, 'StimType', 'laser'); end
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
%       varargin    - 'AlignToStim': whether to align to first stim start
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'StimType': stimulation type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'  - detect automatically
%                       'none'  - no stimulation
%                       'gas'   - hypoxia
%                       'laser' - laser
%                   default == auto
%                   - 'ParseGas': whether to parse gas pulses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ParseLaser': whether to parse laser pulses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpectrogram': whether to plot spectrogram
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotZooms': whether to plot zoomed regions
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotForCorel': whether to plot for CorelDraw
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveTicks': whether to remove all ticks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'StimTableSuffix': Suffix for the stim table
%                   must be a character vector or a string scalar 
%                   default == '_pulses'
%                   - 'RelativeTimeWindow': relative time window
%                   must be a 2-element numeric vector
%                   default == [-10, 20] for gas and [-10, 10] for laser
%                   - 'TimeUnits': output time units
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'min'   - minutes
%                       's'     - seconds
%                       'ms'    - milliseconds
%                       'us'    - microseconds
%                   default == 's'
%                   - 'FigExpansion': expansion factor for figure position
%                   must be a positive scalar or 2-element vector
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
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
%       /home/Matlab/plethR01/plethR01_plot_figureLaserExamples.m

% File History:
% 2019-09-30 Moved from plethRO1_analyze.m
% 2019-10-04 Now computes the spectrogram for the full trace 
%               and then restrict it to the time window of interest
% 

%% Hard-coded constants
S_PER_MIN = 60;

%% Hard-coded parameters
validStimTypes = {'auto', 'none', 'gas', 'laser'};
validTimeUnits = {'', 'min', 's', 'ms', 'us'};
stimTableSuffixGas = '_gas_pulses';
stimTableSuffixLaser = '_laser_pulses';

% TODO: Make optional arguments
spectYLimits = [1, 50];         % Look at 1-50 Hz by default
% relativeTimeWindowNoSwd = [-6.3, -6];
% relativeTimeWindowSwd = [6, 6.3];
suffixNoSwd = '_no_swd';
suffixSwd = '_swd';

%% Default values for optional arguments
alignToStimDefault = false;
stimTypeDefault = 'auto';
parseGasDefault = false;
parseLaserDefault = false;
plotSpectrogramDefault = true;
plotZoomsDefault = false;
plotForCorelDefault = false;
removeTicksDefault = false;
stimTableSuffixDefault = '';    % set later

% TODO: Make optional arguments
channelNamesOriginal = {};
channelNamesToPlot = {};

relativeTimeWindowDefault = [];     % set later
timeUnitsDefault = '';          % set later
figExpansionDefault = [];       % no figure expansion by default
figNameDefault = '';            % don't save figure by default
figTypesDefault = {'png', 'epsc'};

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
addParameter(iP, 'AlignToStim', alignToStimDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'StimTypes', stimTypeDefault, ...
    @(x) any(validatestring(x, validStimTypes)));
addParameter(iP, 'ParseGas', parseGasDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ParseLaser', parseLaserDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSpectrogram', plotSpectrogramDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotZooms', plotZoomsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotForCorel', plotForCorelDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveTicks', removeTicksDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'StimTableSuffix', stimTableSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RelativeTimeWindow', relativeTimeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));
addParameter(iP, 'FigExpansion', figExpansionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, spike2Path, varargin{:});
alignToStim = iP.Results.AlignToStim;
stimType = validatestring(iP.Results.StimTypes, validStimTypes);
parseGas = iP.Results.ParseGas;
parseLaser = iP.Results.ParseLaser;
plotSpectrogram = iP.Results.PlotSpectrogram;
plotZooms = iP.Results.PlotZooms;
plotForCorel = iP.Results.PlotForCorel;
removeTicks = iP.Results.RemoveTicks;
stimTableSuffix = iP.Results.StimTableSuffix;
relativeTimeWindow = iP.Results.RelativeTimeWindow;
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);
figExpansion = iP.Results.FigExpansion;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_traces() function
% otherArguments = iP.Unmatched;

%% Preparation
% Force the Spike2 path to end with '.mat'
spike2Path = force_string_end(spike2Path, '.mat');

% Detect stimulation type
if strcmpi(stimType, 'auto')
    % If not aligning, stim type is 'none'
    if ~alignToStim
        stimType = 'none';
    else
        % Try finding the corresponding gas stim table
        [~, stimTablePath] = ...
            find_matching_files(spike2Path, 'Suffix', stimTableSuffixGas, ...
                                'Extension', 'csv');

        % Otherwise, try finding the corresponding laser stim table
        if ~isempty(stimTablePath)
            stimType = 'gas';        
        else
            [~, stimTablePath] = ...
                find_matching_files(spike2Path, 'Suffix', stimTableSuffixLaser, ...
                                    'Extension', 'csv');
                                
            % If not found, don't align to stimulus
            if ~isempty(stimTablePath)
                stimType = 'laser';
            else
                alignToStim = false;
                stimType = 'none';
            end
        end
    end
elseif strcmpi(stimType, 'gas')
    [~, stimTablePath] = ...
        find_matching_files(spike2Path, 'Suffix', stimTableSuffixGas, ...
                            'Extension', 'csv');
    if isempty(stimTablePath)
        parseGas = true;
    end
elseif strcmpi(stimType, 'laser')
    [~, stimTablePath] = ...
        find_matching_files(spike2Path, 'Suffix', stimTableSuffixLaser, ...
                            'Extension', 'csv');
    if isempty(stimTablePath)
        parseLaser = true;
    end
end

% Decide on stimulation table suffix
if isempty(stimTableSuffix)
    switch stimType
        case 'gas'
            stimTableSuffix = stimTableSuffixGas;
        case 'laser'
            stimTableSuffix = stimTableSuffixLaser;
        otherwise
            stimTableSuffix = '';
    end
end

% Decide on original channel names
if isempty(channelNamesOriginal)
    switch stimType
        case 'gas'
            channelNamesOriginal = {'O2'; 'Pleth 2'; 'WIC#2'; 'WIC#1'};
        case 'laser'
            channelNamesOriginal = {'Sound'; 'Pleth 2'; 'WIC#2'; 'WIC#1'};
        case 'none'
            channelNamesOriginal = {'Pleth 2'; 'WIC#2'; 'WIC#1'};
        otherwise
    end
end

% Decide on channel names to plot
if isempty(channelNamesToPlot)
    switch stimType
        case 'gas'
            channelNamesToPlot = {'O2'; 'Pleth'; 'EEG'; 'EMG'};
        case 'laser'
            channelNamesToPlot = {'Laser'; 'Pleth'; 'EEG'; 'EMG'};
        otherwise
            channelNamesToPlot = {'Pleth'; 'EEG'; 'EMG'};
    end
end

% Decide on relative time window
if isempty(relativeTimeWindow)
    switch stimType
        case 'gas'
            relativeTimeWindow = [-10, 20];
        case 'laser'
            relativeTimeWindow = [-10, 10];
        otherwise
            relativeTimeWindow = [];
    end
end

% Decide on time units
if isempty(timeUnits)
    switch stimType
        case 'gas'
            timeUnits = 'min';
        case 'laser'
            timeUnits = 's';
        otherwise
            timeUnits = 's';
    end
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
if isempty(figExpansion)
    if plotForCorel
        figExpansion = [1, 0.4];
    else
        figExpansion = [1, 1];
    end
end

% Decide on a figure name
if isempty(figName)
    figPathBase = extract_fileparts(spike2Path, 'pathbase');

    if ~isempty(relativeTimeWindow)
        figPathBase = strcat(figPathBase, '_rel', ...
                            num2str(relativeTimeWindow(1)), 'to', ...
                            num2str(relativeTimeWindow(2)));
    end
else
    figPathBase = extract_fileparts(figName, 'pathbase');
end

% Create other figure names
figPathBaseNoSwd = strcat(figPathBase, suffixNoSwd);
figPathBaseSwd = strcat(figPathBase, suffixSwd);

%% Load trace data
% Load data and parse gas trace if necessary
data = parse_spike2_mat(spike2Path, 'ChannelNames', channelNamesOriginal, ...
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

%% Compute the spectrogram from the EEG trace
% Look for an EEG trace
eegChannelNum = find_in_strings('EEG', channelNamesToPlot, 'MaxNum', 1);

% Only compute spectrogram if EEG trace is found
if plotSpectrogram && ~isempty(eegChannelNum)
    % Update channel names to plot
    channelNamesToPlot = [channelNamesToPlot; {'Spect'}];

    % Get the entire EEG trace
    eegChannelValues = channelValues(:, eegChannelNum);

    % Compute the spectrogram
    [spectData, freqHz, timeInstantsSeconds] = ...
        compute_spectrogram(eegChannelValues, siSeconds);

    % Convert times to minutes
    timeInstants = timeInstantsSeconds * unitsPerSecond;
else
    plotSpectrogram = false;
end

%% Restrict or align trace data if requested
if alignToStim
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
    timeWindow = centerTime + relativeTimeWindow;

    % Find the time window endpoints
    endPoints = find_window_endpoints(timeWindow, timeVec);

    % Restrict time vector and shift it to be 
    %   relative to the first stim start time
    timeVecToPlot = extract_subvectors(timeVec, 'EndPoints', endPoints);
    timeVecToPlot = timeVecToPlot - centerTime;

    % Restrict data to the time window endpoints
    channelValuesToPlot = ...
        extract_subvectors(channelValues, 'EndPoints', endPoints);

    % Restrict spectrogram to the time window endpoints
    if plotSpectrogram
        % Find the time window endpoints in timeInstants
        %   making sure that all parts of the time window is included
        endPointsInstants = find_window_endpoints(timeWindow, timeInstants, ...
                                            'BoundaryMode', 'inclusive');

        % Restrict timeInstants to the time window endpoints and shift it to be 
        %   relative to the first stim start time
        timeInstants = ...
            extract_subvectors(timeInstants, 'EndPoints', endPointsInstants);
        timeInstantsToPlot = timeInstants - centerTime;

        % Restrict spectrogram data to the time window endpoints
        spectDataToPlot = transpose(extract_subvectors(transpose(spectData), ...
                                            'EndPoints', endPointsInstants));
    end
else
    timeVecToPlot = timeVec;
    channelValuesToPlot = channelValues;
    if plotSpectrogram
        timeInstantsToPlot = timeInstants;
        spectDataToPlot = spectData;
    end
end

%% Plot traces
% Count the number of subplots
nSubPlots = numel(channelNamesToPlot);

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
        imagesc(ax(iPlot), timeInstantsToPlot, freqHz, abs(spectDataToPlot));

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
    if alignToStim
        plot_vertical_line(0, 'AxesHandle', ax(iPlot), 'Color', 'k', ...
                            'LineWidth', 1);
    end
end

% Link the x axes
linkaxes(ax, 'x');

% Set x axes limits to the relative window
if alignToStim && ~isempty(relativeTimeWindow)
    xlim(relativeTimeWindow);
end

% Update figure for Corel Draw
if plotForCorel
    fig = update_figure_for_corel(fig, 'RemoveTicks', removeTicks);
elseif removeTicks
    % TODO
end

% Save the figure in different zooms
if plotZooms
    % Shrink figure width
    set_figure_properties('FigHandle', fig, 'FigExpansion', [0.5, 1], ...
                            'ExpandFromDefault', false);

    % Zoom to first region and save
    xlim(relativeTimeWindowNoSwd);
    save_all_figtypes(fig, figPathBaseNoSwd, figTypes);

    % Zoom to second region and save
    xlim(relativeTimeWindowSwd);
    save_all_figtypes(fig, figPathBaseSwd, figTypes);

    % Expand figure width
    set_figure_properties('FigHandle', fig, 'FigExpansion', [2, 1], ...
                            'ExpandFromDefault', false);

    % Plot shades
    for iPlot = 1:numel(channelNamesToPlot)
        % Plot vertical shades
        subplot(ax(iPlot))
        plot_window_boundaries([relativeTimeWindowNoSwd, relativeTimeWindowSwd], ...
                                'BoundaryType', 'verticalShade');
    end

    % Zoom back out
    xlim(relativeTimeWindow);
end

% Save figure
save_all_figtypes(fig, figPathBase, figTypes);    

%% Output results
handles.fig = fig;
handles.ax = ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spectData, freqHz, timeInstantsSeconds] = ...
                compute_spectrogram(eegValues, siSeconds)
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
%   Note: time instants are the midpoints of each time window
[spectData, freqHz, timeInstantsSeconds] = ...
    spectrogram(eegValues, binWidthSamples, ...
                overlapSamples, [], samplingFreqHz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

spike2PathBase = extract_fileparts(spike2Path, 'pathbase');
stimTablePath = [spike2PathBase, stimTableSuffix, '.csv'];

% Extract the channel names and values
channelNames = data.channelNames;

% Create figure paths
[figPathBase, figPathBaseNoSwd, figPathBaseSwd] = ...
    argfun(@(x) fullfile(outFolder, x), figBase, figBaseNoSwd, figBaseSwd);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%