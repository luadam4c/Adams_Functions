function [data, siUs, tVec, siPlot] = plot_traces_abf (fileName, varargin)
%% Takes an abf file and plots all traces
% Usage: [data, siUs, tVec, siPlot] = plot_traces_abf (fileName, varargin)
% Explanation:
%       TODO: Uses abf2load, and if not found, uses abfload
% Outputs:
%       data        - full data
%       siUs        - sampling interval in microseconds
%       tVec        - a time vector that can be used to plot things, 
%                       units are in timeUnits (see below for default)
%       siPlot      - sampling interval used for plotting
% Arguments:
%       fileName    - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       varargin    - 'ExpMode': experiment mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'EEG'   - EEG data; x axis in seconds; y-axis in uV
%                       'patch' - patch data; x axis in ms; y-axis in mV
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
%                   default == 's' for EEG data and 'ms' for patch data
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
%                   - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   TODO:
%                   - 'Data': full data
%                   must be a TODO
%                   default == what the file provides
%                   - 'SiUs': sampling interval in microseconds
%                   must be a TODO
%                   default == what the file provides
%
% Requires:
%       cd/check_dir.m
%       cd/parse_abf.m
%       cd/plot_traces.m
%
% Used by:
%       cd/plot_all_abfs_dir.m
%       /media/shareX/share/Adam/Sample_files_from_Katie/test_sweeps.m
%
% File history: 
% 2016-09-22 - adapted from plot_traces_abf_EEG
% 2017-02-13 - BT - added labelling detection between current and voltage
% 2017-03-01 - BT - moved duplicate plotting code into separate helper function
% 2017-03-15 - BT - adapted labelling detection for conductance and differing recmodes
% 2017-04-11 - Fixed the case when range == 0 in plot_traces_abf_helper
% 2017-04-11 - Added outFolder as an optional argument
% 2017-04-13 - Added data, siUs as optional arguments
% 2017-06-16 - Now uses identify_channels.m
% 2017-06-16 - Changed channelLabels to include units
% 2018-01-24 - Added isdeployed
% 2018-07-24 - Now uses a try catch statement
% 2018-09-17 - Added the input parser
% 2018-09-17 - Moved code to parse_abf.m
% 2018-09-18 - Renamed plot_traces_abf_helper -> plot_traces
%               and moved to its own function
% TODO: Improve code legibility with usage of dataReordered instead of data
%

%% Hard-coded parameters
validExpModes = {'EEG', 'patch', ''};
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Undefined'};

%% Default values for optional arguments
expModeDefault = 'patch';       % assume traces are patching data by default
individuallyDefault = false;    % plot all sweeps together by default
outFolderDefault = '';          % set later
timeUnitsDefault = '';          % set later
timeStartDefault = [];          % set later
timeEndDefault = [];            % set later
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later
verboseDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'fileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ExpMode', expModeDefault, ...
    @(x) any(validatestring(x, validExpModes)));
addParameter(iP, 'Individually', individuallyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) isempty(x) || isnumeric(x) && isscalar(x) && x >= 0);
addParameter(iP, 'TimeEnd', timeEndDefault, ...
    @(x) isempty(x) || isnumeric(x) && isscalar(x) && x >= 0);
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || iscellstr(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || iscellstr(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, fileName, varargin{:});
expMode = validatestring(iP.Results.ExpMode, validExpModes);
individually = iP.Results.Individually;
outFolder = iP.Results.OutFolder;
timeUnits = iP.Results.TimeUnits;
timeStart = iP.Results.TimeStart;
timeEnd = iP.Results.TimeEnd;
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;
verbose = iP.Results.Verbose;

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

% Set (some) dependent argument defaults
[fileDir, fileBase, ~] = fileparts(fileName);
if isempty(outFolder)
    outFolder = fullfile(fileDir, strcat(fileBase, '_traces'));
end
if verbose
    fprintf('Outfolder is %s ...\n', outFolder);
end

%% Check if needed output directories exist
check_dir(outFolder, 'Verbose', verbose);

%% Load data
% Load and parse the abf file
[abfParams, data, tVec, ~, ~, ~, dataReordered] = ...
    parse_abf(fileName, 'Verbose', false, ...
              'ExpMode', expMode, 'TimeUnits', timeUnits, ...
              'ChannelTypes', channelTypes, ...
              'ChannelUnits', channelUnits, ...
              'ChannelLabels', channelLabels);

% Extract the parsed parameters
expMode = abfParams.expMode;
channelLabels = abfParams.channelLabels;
nDimensions = abfParams.nDimensions;
nChannels = abfParams.nChannels;
nSweeps = abfParams.nSweeps;
siUs = abfParams.siUs;
siPlot = abfParams.siPlot;
timeUnits = abfParams.timeUnits;

%% Prepare for plotting
% Set default x-axis limits
if isempty(timeStart)
    timeStart = 0;
end
if isempty(timeEnd)
    timeEnd = tVec(end);
end
if verbose
    fprintf('Interval to show = [%g %g]\n', timeStart, timeEnd);
end
xlimits = [timeStart, timeEnd];

% Set up labels for each trace based on experiment mode
switch expMode
case 'EEG'
    traceLabels = cell(1, nChannels);
    for iChannel = 1:nChannels
        traceLabels{iChannel} = ['Channel #', num2str(iChannel)];
    end
case 'patch'
    traceLabels = cell(1, nSweeps);
    for iSwp = 1:nSweeps
        traceLabels{iSwp} = ['Sweep #', num2str(iSwp)];
    end
otherwise
    error('Expmode unrecognized!');
end

% Set up the time axis label
xLabel = ['Time (', timeUnits, ')'];

%% Do the plotting
% Construct a file identifier
fileIdentifier = sprintf('%s_%g%s_%g%s', fileBase, xlimits(1), timeUnits, ...
                                        xlimits(2), timeUnits);

% Plot data
if ~individually && strcmpi(expMode, 'EEG')
    % Print message
    if verbose
        fprintf('Plotting all channels ...\n');
    end

    % Decide on figure name and title
    vecAll = data;
    yLabel = channelLabels{1};
    figTitle = sprintf('All channels for %s', fileIdentifier);
    figName = fullfile(outFolder, sprintf('%s_all.png', fileIdentifier));
    figNum = 1;

    % Do the plotting
    h = plot_traces(tVec, vecAll, xlimits, xLabel, yLabel, ...
                    traceLabels, figTitle, figName, figNum);
elseif ~individually && strcmpi(expMode, 'patch') || ...
        individually && strcmpi(expMode, 'EEG')
    % Loop through all channels
    parfor iChannel = 1:nChannels
        % Print message
        if verbose
            fprintf('Plotting all sweeps of Channel #%d ...\n', iChannel);
        end

        % Decide on figure name and title
        if nDimensions == 2
            vecAll = data(:, iChannel);
        else
            vecAll = squeeze(data(:, iChannel, :));
        end
        yLabel = channelLabels{iChannel};
        figTitle = sprintf('Data for Channel #%d for %s', ...
                            iChannel, fileIdentifier);
        figName = fullfile(outFolder, sprintf('%s_Channel%d_all.png', ...
                                            fileIdentifier, iChannel));
        figNum = 100 * iChannel;

        % Do the plotting
        h = plot_traces(tVec, vecAll, xlimits, xLabel, yLabel, ...
                        traceLabels, figTitle, figName, figNum);
    end
elseif individually && strcmpi(expMode, 'patch')
    % Plot each channel and sweep individually
    for iChannel = 1:nChannels
        parfor iSwp = 1:nSweeps
            % Print message
            if verbose
                fprintf('Plotting Channel #%d, Sweep #%d ...\n', ...
                            iChannel, iSwp);
            end

            % Decide on figure name and title
            if nDimensions == 2
                vecAll = data(:, iChannel);
            elseif nDimensions == 3
                vecAll = data(:, iChannel, iSwp);
            end
            yLabel = channelLabels{iChannel};
            figTitle = sprintf('Data for Channel #%d, Sweep #%d for %s', ...
                                iChannel, iSwp, fileIdentifier);
            figName = fullfile(outFolder, ...
                                sprintf('%s_Channel%d_Sweep%d.png', ...
                                        fileIdentifier, iChannel, iSwp));
            figNum = 100 * iChannel + iSwp;

            % Do the plotting
            h = plot_traces(tVec, vecAll, xlimits, xLabel, yLabel, ...
                            traceLabels, figTitle, figName, figNum);
        end
    end
else
    error('Not Implemented Yet!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
OLD CODE:

%    minimum = min(min(vvecAll));
%    maximum = max(max(vvecAll));
%    range = maximum - minimum;
%    h = figure(1);
%    set(h, 'Visible', 'Off');
%    clf(h);
%    for j = 1:nChannels
%        plot(tVec, vvecAll(:, j));    hold on;
%    end
%    axis([timeStart timeEnd minimum-0.2*range maximum+0.2*range]);
%    xlabel(['Time (', timeUnits, ')']);
%
%        minimum = min(min(vecAll));
%        maximum = max(max(vecAll));
%        range = maximum - minimum;
%        h = figure(100*j);
%        set(h, 'Visible', 'Off');
%        clf(h);
%        for k = 1:nSweeps
%            plot(tVec, vecAll(:, k));    hold on;
%        end
%        axis([timeStart timeEnd minimum-0.2*range maximum+0.2*range]);
%        xlabel(['Time (', timeUnits, ')']);
%
%        minimum = min(vvec);
%        maximum = max(vvec);
%        range = maximum - minimum;
%        h = figure(1+j);
%        set(h, 'Visible', 'Off');
%        clf(h);
%        plot(tVec, vvec, 'k');    hold on;
%        axis([timeStart timeEnd minimum-0.2*range maximum+0.2*range]);
%        xlabel(['Time (', timeUnits, ')']);
%
%            minimum = min(vec);
%            maximum = max(vec);
%            range = maximum - minimum;
%            h = figure(100*j + k);
%            set(h, 'Visible', 'Off');
%            clf(h);
%            plot(tVec, vec, 'k');    hold on;
%            axis([timeStart timeEnd minimum-0.2*range maximum+0.2*range]);
%            xlabel(['Time (', timeUnits, ')']);
%         % y-axis labels: channelLabels{1} is Voltage, channelLabels{2} is Current, channelLabels{3} is Conductance
%        if j == ind_current    % if the current channel is "Current"
%            ylabel(sprintf('%s (%s)', channelLabels{2}, channelUnits{2}));        
%        else            %%% what happens if we have conductance channels too?
%            ylabel(sprintf('%s (%s)', channelLabels{1}, channelUnits{1}));
%        end
%if strcmp(expMode, 'EEG')
%    if nDimensions == 1
%        channelUnits = 'mV';
%    elseif nDimensions == 2
%        channelUnits = {'uV', 'uV'};
%    elseif nDimensions == 3
%        channelUnits = {'uV', 'uV', 'uV'};
%    end
%    channelLabels = 'EEG amplitude';
%elseif strcmp(expMode, 'patch')
%    [channelUnits, channelLabels] = correct_patch_labels(nDimensions, data);
%end
%


    if nDimensions == 2        % Usually EEG         %TODO: Not always true, should be fixed
        timeUnits = 's';
    elseif nDimensions == 3    % Usually Patch clamp    %TODO: Not always true, should be fixed
        timeUnits = 'ms';
    end
    if nDimensions == 2        % Usually EEG
        channelUnits = 'uV';
    elseif nDimensions == 3    % Usually Patch clamp
        channelUnits = {'mV', 'pA', 'nS'};
    end
    if nDimensions == 2        % Usually EEG
        channelLabels = 'EEG amplitude';
    elseif nDimensions == 3    % Usually Patch clamp
        channelLabels = {'Voltage', 'Current', 'Conductance'};
    end

%    recmode        - (opt)    1: (mV)        %%% BT: recmode never explicitly implemented, documentation for possible channel combinations
%                2: (uV, uV)
%                3: (mV, pA)
%                4: (pA, mV)
%                5: (uV, uV, uV)
%                6: (mV, pA, nS)
%                7: (mV, mV, pA)
%                8: (pA, mV, nS)
%                sampling interval is assumed to be in microseconds
%                default == TODO

function [channelUnits, channelLabels] = correct_patch_labels(data)
%% Finds correct units and labels of channels if data is patch clamp
% Usage: [channelUnits, channelLabels] = correct_patch_labels(data)
%    data        - raw data vector
% Used by:
%        /media/shareX/brianT/ABF_plotting/plot_traces_abf_bt.m

nChannels = size(data, 2);
channelUnits = cell(1, nChannels);
channelLabels = cell(1, nChannels);
ranges = zeros(1, nChannels);    % stores maximum absolute range for each channel
peak = zeros(1, nChannels);    % maximum for each channel
avgs = zeros(1, nChannels);    % average for each channel
for j = 1:nChannels        % for each channel
    vecAll = squeeze(data(:, j, :));                % all traces in this channel
    ranges(j) = abs(max(max(vecAll)) - min(min(vecAll)));    % maximum range of all traces in this channel
    peak(j) = max(max(vecAll));
end
avgs = mean(mean(data, 1), 3);
if nChannels == 2
    [~, ind_current] = max(ranges);
    ranges(ind_current) = -Inf;
    [~, ind_voltage] = max(ranges);
    channelUnits{ind_voltage} = 'mV';
    channelUnits{ind_current} = 'pA';
    channelLabels{ind_voltage} = 'Voltage';
    channelLabels{ind_current} = 'Current';
elseif nChannels == 3
    if abs(range(1) - range(2)) < 5    % if ranges of first two channels are similar enough, assume both are voltage
        channelUnits = {'mV', 'mV', 'pA'};         % default: channelUnits = {'mV', 'pA', 'nS'};
        channelLabels = {'Voltage', 'Voltage', 'Current'};    % default: channelLabels = {'Voltage', 'Current', 'Conductance'};
    else                
        [~, ind_conductance] = max(avgs);    % greatest average is presumably "Conductance" (generally positive)
        ranges(ind_conductance) = -Inf;        % all ranges are positive, remove conductance from current detection
        [~, ind_current] = max(ranges);        % the channel with largest range is presumably "Current"
        ranges(ind_current) = -Inf;        % remove current from ranges, leaving only "Voltage"
        [~, ind_voltage] = max(ranges);
        channelUnits{ind_voltage} = 'mV';        % Reorder channelUnits and channelLabels to new indices
        channelUnits{ind_current} = 'pA';
        channelLabels{ind_voltage} = 'Voltage';
        channelLabels{ind_current} = 'Current';
        channelLabels{ind_conductance} = 'Conductance';
        if peak(ind_conductance) > 2000    % Conductance usually exceeds 2000 when pico    
            channelUnits{ind_conductance} = 'pS';
        else
            channelUnits{ind_conductance} = 'nS';
        end
    end
end

        [channelUnits, ~] = correct_patch_labels(data);
        [~, channelLabels] = correct_patch_labels(data);

ylabel(sprintf('%s (%s)', channelLabels, channelUnits));
ylabel(sprintf('%s (%s)', channelLabels{j}, channelUnits{j}));

if nargin < 1
    error('No fileName specified!');
elseif ~ischar(fileName)
    error('Filename must be a char array in single quotes!');
elseif nargin < 2
    error('No expMode specified!');
elseif ~strcmp(expMode, 'EEG') && ~strcmp(expMode, 'patch')
    error('Expmode must be ''EEG'' or ''patch''!');
elseif nargin >= 5 && ~isdir(outFolder)
    error('outFolder must be a directory!');
elseif nargin >= 8 && plotMode ~= 1 && plotMode ~= 2
    error('The plot mode %d is not currently supported!', plotMode);
elseif nargin >= 9 && ~ischar(timeUnits)
    error('timeUnits must be a char array!');
end

% Set default plot mode
if nargin < 8 || isempty(plotMode)
    plotMode = 1;
end
fprintf('Using plot mode == %d ...\n', plotMode);

if nargin < 3
end

if nargin < 11 || isempty(channelUnits)
    if strcmp(expMode, 'EEG')
        if size(data, 2) == 1
            channelUnits = 'mV';
        else
            channelUnits = 'uV';
        end
    elseif strcmp(expMode, 'patch')
        [~, channelUnits, ~] = identify_channels(data);
    end
elseif nDimensions == 2 && ~ischar(channelUnits) ...
    || nDimensions == 3 && ~iscellstr(channelUnits)
    error('channelUnits must be a char array for 2-data data and a cell array of char arrays for 3-data data!');
end
if nargin < 11 || isempty(channelLabels)
    if strcmp(expMode, 'EEG')
        channelLabels = 'EEG amplitude (mV)';
    elseif strcmp(expMode, 'patch')
        [~, ~, channelLabels] = identify_channels(data);
    end
elseif nDimensions == 2 && ~ischar(channelLabels) ...
    || nDimensions == 3 && ~iscellstr(channelLabels)
    error('channelLabels must be a char array for 2-data data and a cell array of char arrays for 3-data data!');
end

%   TODO: update below:
%       channelUnits    - (opt) units for each channel, must be a char array for 2-data data 
%                   and a cell array of char arrays for 3-data data
%                   default == 'uV' for 2-data data and {'mV', 'pA', 'nS'} for 3-data data)
%       channelLabels   - (opt) labels for each channel, must be a char array for 2-data data
%                   and a cell array of char arrays for 3-data data
%                   default == 'EEG amplitude' for 2-data data and {'Voltage', 'Current', 'Conductance'} for 3-data data

function [data, siUs, tVec] = plot_traces_abf (fileName, expMode, data, siUs, outFolder, timeStart, timeEnd, plotMode, timeUnits, channelUnits, channelLabels)

if strcmp(timeUnits, 'ms')
    % Use a sampling interval in ms
    siPlot = siUs / MS_PER_S;
elseif strcmp(timeUnits, 's')
    % Use a sampling interval in seconds
    siPlot = siUs / US_PER_S;
end
fprintf('Sampling interval = %d %s\n', siPlot, timeUnits);

if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end

if isempty(timeUnits)
    switch expMode
    case 'EEG'
        timeUnits = 's';
    case 'patch'
        timeUnits = 'ms';
    otherwise
        error('ExpMode unrecognize!');
    end
end

axis([timeStart, timeEnd, minY - 0.2 * rangeY, maxY + 0.2 * rangeY]);

% Used by:
%        /media/shareX/brianT/ABF_plotting/plot_traces_abf_bt.m

vvecAll = data;

%    timeStart  - time interval start
%    timeEnd    - time interval end
xlim([timeStart, timeEnd]);

% Set up labels for each trace
if nDimensions == 2        % Usually EEG
    traceLabels = cell(1, nChannels);
    for j = 1:nChannels
        traceLabels{j} = ['Channel #', num2str(j)];
    end
elseif nDimensions == 3    % Usually Patch clamp
    traceLabels = cell(1, nSweeps);
    for k = 1:nSweeps
        traceLabels{k} = ['Sweep #', num2str(k)];
    end
end

plot(tVec, dataVec(:, iTrace), ...
    'DisplayName', ['Sweep #', num2str(iTrace)]);

legend(traceLabels);

fileIdentifier = sprintf('%s_%.1f_%.1f', fileBase, xlimits(1), xlimits(2));
title(sprintf('Data for all channels between %.1f %s and %.1f %s', timeStart, timeUnits, timeEnd, timeUnits));

%    timeUnits  - units of time
%    nTraces    - number of channels if plotMode == 1, 1 if plotMode == 2

title(sprintf('Data for Channel #%d between %.1f %s and %.1f %s', ...
                iChannel, timeStart, timeUnits, timeEnd, timeUnits));

% Plot raw data (each channel and/or sweep individually)
if nDimensions == 2 && individually        % could be EEG or patch clamp
    for iChannel = 1:nChannels
        if verbose
            fprintf('Plotting channel #%d ...\n', iChannel);
        end
        figName = fullfile(outFolder, ...
                    sprintf('%s_Channel%d.png', ...
                    fileIdentifier, iChannel))
        vvec = data(:, iChannel);
        figNum = 1+iChannel;
        h = plot_traces(figNum, tVec, vvec, ...
                                    xlimits, xLabel, yLabel, ...
                                    traceLabels, figTitle, figName);
        title(sprintf('Data for %s between %.1f %s and %.1f %s', traceLabels{iChannel}, timeStart, timeUnits, timeEnd, timeUnits));
        ylabel(channelLabels);
        saveas(h, figName, 'png');
        hold off;
        close(h);
    end
elseif nDimensions == 3 && individually    % usually Patch clamp        %%% Need to fix this part too
    for iChannel = 1:nChannels
        for iSwp = 1:nSweeps
            if verbose
                fprintf('Plotting channel #%d and sweep #%d ...\n', iChannel, iSwp);
            end
            figName = fullfile(outFolder, ...
                                sprintf('%s_Channel%d_Sweep%d.png', ...
                                        fileIdentifier, iChannel, iSwp));
            vec = data(:, iChannel, iSwp);
            figNum = 100*iChannel+iSwp;
            h = plot_traces(figNum, tVec, vec, ...
                                       xlimits, xLabel, yLabel, ...
                                       traceLabels, figTitle, figName);
            title(sprintf('Data for %s between %.1f %s and %.1f %s', ...
                    traceLabels{iSwp}, timeStart, timeUnits, timeEnd, timeUnits));
            ylabel(channelLabels);
            saveas(h, figName, 'png');
            hold off;
            close(h);
        end
    end
end

%}
