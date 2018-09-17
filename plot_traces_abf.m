function [data, siUs, tVec] = plot_traces_abf (fileName, varargin)
%% Takes an abf file and plots all traces
% Usage: [data, siUs, tVec] = plot_traces_abf (fileName, varargin)
% Explanation:
%       TODO: Uses abf2load, and if not found, uses abfload
% Outputs:
%       data        - full data
%       siUs        - sampling interval in microseconds
%       tVec        - a time vector that can be used to plot things, 
%                       units are in timeUnits (see below for default)
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
%                   default == 's' for 2-data data and 'ms' for 3-data data
%                   - 'TimeStart': the start of the time interval of interest 
%                                   (in units set by TimeUnits)
%                   must be a numeric nonnegative scalar
%                   default == 0
%                   - 'TimeEnd': the end of the time interval of interest 
%                                   (in units set by TimeUnits)
%                   must be a numeric nonnegative scalar
%                   default == tVec(end)
%                   TODO:
%                   - 'Data': full data
%                   must be a TODO
%                   default == what the file provides
%                   - 'SiUs': sampling interval in microseconds
%                   must be a TODO
%                   default == what the file provides
%
% Requires:
%       cd/parse_abf.m
%
% Used by:
%       /media/shareX/share/Adam/Sample_files_from_Katie/test_sweeps.m
%
% File history: 
% 2016-09-22 - adapted from plot_traces_abf_EEG
% 2017-02-13 - BT - added labelling detection between current and voltage
% 2017-03-01 - BT - moved duplicate plotting code into separate helper function
% 2017-03-15 - BT - adapted labelling detection for conductance and differing recmodes
% 2017-04-11 - Fixed the case when range == 0 in plot_traces_abf_bt_helper
% 2017-04-11 - Added outFolder as an optional argument
% 2017-04-13 - Added data, siUs as optional arguments
% 2017-06-16 - Now uses identify_channels.m
% 2017-06-16 - Changed channelLabels to include units
% 2018-01-24 - Added isdeployed
% 2018-07-24 - Now uses a try catch statement
% 2018-09-17 - Added the input parser
% 2018-09-17 - Moved code to parse_abf.m
%

%% Hard-coded parameters
validExpModes = {'EEG', 'patch'};

%% Default values for optional arguments
expModeDefault = 'patch';       % assume traces are patching data by default
individuallyDefault = false;    % plot all sweeps together by default
outFolderDefault = '';          % set later
timeUnitsDefault = '';          % set later
timeStartDefault = [];          % set later
timeEndDefault = [];            % set later

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

% Read from the Input Parser
parse(iP, fileName, varargin{:});
expMode = validatestring(iP.Results.ExpMode, validExpModes);
individually = iP.Results.Individually;
outFolder = iP.Results.OutFolder;
timeUnits = iP.Results.TimeUnits;
timeStart = iP.Results.TimeStart;
timeEnd = iP.Results.TimeEnd;

% Set (some) dependent argument defaults
if isempty(outFolder)
    [fileDir, fileBase, ~] = fileparts(fileName);
    outFolder = fullfile(fileDir, strcat(fileBase, '_traces'));
end
fprintf('Outfolder is %s ...\n', outFolder);
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

%% Check if needed output directories exist
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end

%% Load data and prepare for plotting
% Load and parse the abf file
[data, siUs, abfParams] = parse_abf(fileName);

% Extract the parsed parameters
channelUnits = abfParams.channelUnits;
channelLabels = abfParams.channelLabels;
ndim = abfParams.ndim;
nSamples = abfParams.nSamples;
nChannels = abfParams.nChannels;
nSweeps = abfParams.nSweeps;

% Get the sampling interval for plotting
if strcmp(timeUnits, 'ms')
    % Use a sampling interval in ms
    si = abfParams.siMs;
elseif strcmp(timeUnits, 's')
    % Use a sampling interval in seconds
    si = abfParams.siSeconds;
end

% Construct time vector
tVec = si * (1:nSamples)';

% Set default x-axis limits
if isempty(timeStart)
    timeStart = 0;
end
if isempty(timeEnd)
    timeEnd = tVec(end);
end
fprintf('Interval to show = [%g %g]\n', timeStart, timeEnd);

% Set up labels for each trace
if ndim == 2        % Usually EEG
    traceLabels = cell(1, nChannels);
    for j = 1:nChannels
        traceLabels{j} = ['Channel #', num2str(j)];
    end
elseif ndim == 3    % Usually Patch clamp
    traceLabels = cell(1, nSweeps);
    for k = 1:nSweeps
        traceLabels{k} = ['Sweep #', num2str(k)];
    end
end

%% Do the plotting
% Plot raw data all at once
if ndim == 2 && ~individually       % could be EEG or patch clamp
    if strcmp(expMode, 'EEG')
        fprintf('Plotting all channels ...\n');
        vvecAll = data;
        figureNum = 1;
        h = plot_traces_abf_bt_helper(figureNum, tVec, vvecAll, timeStart, timeEnd, timeUnits, nChannels);
        title(sprintf('Data for all channels between %.1f %s and %.1f %s', timeStart, timeUnits, timeEnd, timeUnits));
        ylabel(channelLabels);
        legend(traceLabels);
        saveas(h, fullfile(outFolder, sprintf('%.1f_%.1f_all.png', timeStart, timeEnd)), 'png');    
        hold off;
    elseif strcmp(expMode, 'patch')
        % Plot sweeps
        for j = 1:nChannels
            fprintf('Plotting channel #%d for all sweeps ...\n', j);
            vecAll = squeeze(data(:, j, :));
            figureNum = 1*j;
            h = plot_traces_abf_bt_helper(figureNum, tVec, vecAll, timeStart, timeEnd, timeUnits, 1);
            title(sprintf('Data for all channels between %.1f %s and %.1f %s', timeStart, timeUnits, timeEnd, timeUnits));
            ylabel(channelLabels{j});
            legend(traceLabels);
            saveas(h, fullfile(outFolder, sprintf('%.1f_%.1f_Channel%d_all.png', timeStart, timeEnd, j)), 'png');
            hold off;
    %        close(h);
        end
    end
%    close(h);
elseif ndim == 3 && ~individually   % usually Patch clamp
    % Plot sweeps
    for j = 1:nChannels
        fprintf('Plotting channel #%d for all sweeps ...\n', j);
        vecAll = squeeze(data(:, j, :));
        figureNum = 100*j;
        h = plot_traces_abf_bt_helper(figureNum, tVec, vecAll, timeStart, timeEnd, timeUnits, nSweeps);
        title(sprintf('Data for all sweeps between %.1f %s and %.1f %s', timeStart, timeUnits, timeEnd, timeUnits));
        ylabel(channelLabels{j});
        legend(traceLabels);
        saveas(h, fullfile(outFolder, sprintf('%.1f_%.1f_Channel%d_all.png', timeStart, timeEnd, j)), 'png');
        hold off;
%        close(h);
    end
end

% Plot raw data (each channel and/or sweep individually)
if ndim == 2 && individually        % could be EEG or patch clamp
    if strcmp(expMode, 'EEG')
        for j = 1:nChannels
            fprintf('Plotting channel #%d ...\n', j);
            vvec = data(:, j);
            figureNum = 1+j;
            h = plot_traces_abf_bt_helper(figureNum, tVec, vvec, timeStart, timeEnd, timeUnits, 1);
            title(sprintf('Data for %s between %.1f %s and %.1f %s', traceLabels{j}, timeStart, timeUnits, timeEnd, timeUnits));
            ylabel(channelLabels);
            saveas(h, fullfile(outFolder, sprintf('%.1f_%.1f_Channel%d.png', timeStart, timeEnd, j)), 'png');
            hold off;
            close(h);
        end
    elseif strcmp(expMode, 'patch')
        for j = 1:nChannels
            fprintf('Plotting channel #%d ...\n', j);
            vvec = data(:, j);
            figureNum = 1+j;
            h = plot_traces_abf_bt_helper(figureNum, tVec, vvec, timeStart, timeEnd, timeUnits, 1);
            title(sprintf('Data for %s between %.1f %s and %.1f %s', traceLabels{j}, timeStart, timeUnits, timeEnd, timeUnits));
            ylabel(channelLabels{j});
            saveas(h, fullfile(outFolder, sprintf('%.1f_%.1f_Channel%d.png', timeStart, timeEnd, j)), 'png');
            hold off;
            close(h);
        end
    end
elseif ndim == 3 && individually    % usually Patch clamp        %%% Need to fix this part too
    if strcmp(expMode, 'EEG')
        for j = 1:nChannels
            for k = 1:nSweeps
                fprintf('Plotting channel #%d and sweep #%d ...\n', j, k);
                vec = data(:, j, k);
                figureNum = 100*j+k;
                h = plot_traces_abf_bt_helper(figureNum, tVec, vec, timeStart, timeEnd, timeUnits, 1);
                title(sprintf('Data for %s between %.1f %s and %.1f %s', ...
                        traceLabels{k}, timeStart, timeUnits, timeEnd, timeUnits));
                ylabel(channelLabels);
                saveas(h, fullfile(outFolder, sprintf('%.1f_%.1f_Channel%d_Sweep%d.png', ...
                        timeStart, timeEnd, j, k)), 'png');
                hold off;
                close(h);
            end
        end
    elseif strcmp(expMode, 'patch')
        for j = 1:nChannels
            for k = 1:nSweeps
                fprintf('Plotting channel #%d and sweep #%d ...\n', j, k);
                vec = data(:, j, k);
                figureNum = 100*j+k;
                h = plot_traces_abf_bt_helper(figureNum, tVec, vec, timeStart, timeEnd, timeUnits, 1);
                title(sprintf('Data for %s between %.1f %s and %.1f %s', ...
                    traceLabels{k}, timeStart, timeUnits, timeEnd, timeUnits));
                ylabel(channelLabels{j});
                saveas(h, fullfile(outFolder, sprintf('%.1f_%.1f_Channel%d_Sweep%d.png', ...
                    timeStart, timeEnd, j, k)), 'png');
                hold off;
                close(h);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h] = plot_traces_abf_bt_helper(fig_num, tVec, vvec, timeStart, timeEnd, timeUnits, channels_or_sweep_num)
%% Plots traces, sets axes and xlabel
% Usage: [h] = plot_traces_abf_bt_helper(fig_num, tVec, vvec, timeStart, timeEnd, timeUnits, channels_or_sweep_num)
%    h            - figure
%    fig_num            - figure of number
%    tVec            - time vector for plotting
%    vvec            - compressed data vector
%     timeStart            - time interval start
%    timeEnd            - time interval end
%    timeUnits            - units of time
%    channels_or_sweep_num    - number of channels if plotMode == 1, 1 if plotMode == 2
% Used by:
%        /media/shareX/brianT/ABF_plotting/plot_traces_abf_bt.m

minimum = min(min(vvec));
maximum = max(max(vvec));
range = maximum - minimum;
h = figure(fig_num);
set(h, 'Visible', 'Off');
clf(h);
for k = 1:channels_or_sweep_num
    plot(tVec, vvec(:, k));    hold on;
end
if range ~= 0
    axis([timeStart, timeEnd, minimum - 0.2 * range, maximum + 0.2 * range]);
else
    xlim([timeStart, timeEnd]);
end
xlabel(['Time (', timeUnits, ')']);

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
%    if ndim == 1
%        channelUnits = 'mV';
%    elseif ndim == 2
%        channelUnits = {'uV', 'uV'};
%    elseif ndim == 3
%        channelUnits = {'uV', 'uV', 'uV'};
%    end
%    channelLabels = 'EEG amplitude';
%elseif strcmp(expMode, 'patch')
%    [channelUnits, channelLabels] = correct_patch_labels(ndim, data);
%end
%


    if ndim == 2        % Usually EEG         %TODO: Not always true, should be fixed
        timeUnits = 's';
    elseif ndim == 3    % Usually Patch clamp    %TODO: Not always true, should be fixed
        timeUnits = 'ms';
    end
    if ndim == 2        % Usually EEG
        channelUnits = 'uV';
    elseif ndim == 3    % Usually Patch clamp
        channelUnits = {'mV', 'pA', 'nS'};
    end
    if ndim == 2        % Usually EEG
        channelLabels = 'EEG amplitude';
    elseif ndim == 3    % Usually Patch clamp
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
elseif ndim == 2 && ~ischar(channelUnits) ...
    || ndim == 3 && ~iscellstr(channelUnits)
    error('channelUnits must be a char array for 2-data data and a cell array of char arrays for 3-data data!');
end
if nargin < 11 || isempty(channelLabels)
    if strcmp(expMode, 'EEG')
        channelLabels = 'EEG amplitude (mV)';
    elseif strcmp(expMode, 'patch')
        [~, ~, channelLabels] = identify_channels(data);
    end
elseif ndim == 2 && ~ischar(channelLabels) ...
    || ndim == 3 && ~iscellstr(channelLabels)
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
    si = siUs / MS_PER_S;
elseif strcmp(timeUnits, 's')
    % Use a sampling interval in seconds
    si = siUs / US_PER_S;
end
fprintf('Sampling interval = %d %s\n', si, timeUnits);


%}
