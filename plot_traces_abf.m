function [d, sius, tvec] = plot_traces_abf (filename, expmode, d, sius, outfolder, left, right, plotmode, tUnits, chUnits, chLabels)
%% Takes an abf file and plots all traces
% Usage: [d, sius, tvec] = plot_traces_abf (filename, expmode, d, sius, outfolder, left, right, plotmode, tUnits, chUnits, chLabels)
% Explanation:
%   TODO: Uses abf2load, and if not found, uses abfload
% Outputs:
%       d           - full data
%       sius        - sampling interval in microseconds
%       tvec        - a time vector that can be used to plot things, 
%                       units are in tUnits (see below for default)
% Arguments:
%       filename    - must be either the full address or must be in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%       expmode     -   'EEG'
%                       'patch'
%       d           - (opt) full data
%       sius        - (opt) sampling interval in microseconds
%       outfolder   - (opt) the name of the directory that the plots will be placed
%                   must be a character array
%                   default == a subdirectory named by {filename}_traces
%       left        - (opt) the start of the time interval of interest 
%                   (in seconds for 2-d data and in ms for 3-d data)
%                   default == 0
%       right       - (opt) the end of the time interval of interest (in seconds for 2-d data and in ms for 3-d data)
%                   default == tvec(end)
%       plotmode    - (opt) 1: all traces are to be plotted together;
%                           2: each trace is to be plotted individually
%                   default == 1
%       tUnits     - (opt) units for time, must be a character array
%                   default == 's' for 2-d data and 'ms' for 3-d data
%   TODO: update below:
%       chUnits    - (opt) units for each channel, must be a char array for 2-d data 
%                   and a cell array of char arrays for 3-d data
%                   default == 'uV' for 2-d data and {'mV', 'pA', 'nS'} for 3-d data)
%       chLabels   - (opt) labels for each channel, must be a char array for 2-d data
%                   and a cell array of char arrays for 3-d data
%                   default == 'EEG amplitude' for 2-d data and {'Voltage', 'Current', 'Conductance'} for 3-d data
%
% Requires:
%        cd/construct_abffilename.m
%        /home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%        /home/Matlab/Brians_Functions/identify_channels.m
% Used by:
%        /media/shareX/share/Adam/Sample_files_from_Katie/test_sweeps.m
%
% File history: 
% 2016-09-22 - adapted from plot_traces_abf_EEG
% 2017-02-13 - BT - added labelling detection between current and voltage
% 2017-03-01 - BT - moved duplicate plotting code into separate helper function
% 2017-03-15 - BT - adapted labelling detection for conductance and differing recmodes
% 2017-04-11 - Fixed the case when range == 0 in plot_traces_abf_bt_helper
% 2017-04-11 - Added outfolder as an optional argument
% 2017-04-13 - Added d, sius as optional arguments
% 2017-06-16 - Now uses identify_channels.m
% 2017-06-16 - Changed chLabels to include units
% 2018-01-24 - Added isdeployed

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
    addpath(fullfile(functionsdirectory, '/Brians_Functions/'));    
                                                    % for identify_channels.m
end

%% Check arguments
if nargin < 1
    error('No filename specified!');
elseif ~ischar(filename)
    error('Filename must be a char array in single quotes!');
elseif nargin < 2
    error('No expmode specified!');
elseif ~strcmp(expmode, 'EEG') && ~strcmp(expmode, 'patch')
    error('Expmode must be ''EEG'' or ''patch''!');
elseif nargin >= 5 && ~isdir(outfolder)
    error('outfolder must be a directory!');
elseif nargin >= 8 && plotmode ~= 1 && plotmode ~= 2
    error('The plot mode %d is not currently supported!', plotmode);
elseif nargin >= 9 && ~ischar(tUnits)
    error('tUnits must be a char array!');
end

%% Set defaults
if nargin < 3
    % Load abf file, si is in us
    abffilename_full = construct_abffilename(filename);    % creates full path to abf file robustly
    if exist('abf2load', 'file') == 2
        try
            [d, sius] = abf2load(abffilename_full);
        catch ME
            printf('The file %s cannot be read!\n', abffilename_full);
            rethrow(ME)
            return
        end
    elseif exist('abfload', 'file') == 2
        [d, sius] = abfload(abffilename_full);
    end
end

% Set default output folder to place figures
if nargin < 5 || isempty(outfolder) 
    [filepath, filebase, ~] = fileparts(filename);
    outfolder = fullfile(filepath, strcat(filebase, '_traces'));
end
fprintf('Outfolder is %s ...\n', outfolder);
% Set default plot mode
if nargin < 8 || isempty(plotmode)
    plotmode = 1;
end
fprintf('Using plot mode == %d ...\n', plotmode);

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));    % for abf2load.m or abfload.m

%% Create outfolder if not already exists
if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
    fprintf('New directory made: %s\n\n', outfolder);
end


%% Find data dimensions
ndim = ndims(d);    % number of dimensions in data
if ndim > 3
    error('Cannot parse data with more than 3 dimensions!');
end
fprintf('Number of data dimensions =  %d\n\n', ndim);

%% Units and labels
if nargin < 9 || isempty(tUnits)
    if strcmp(expmode, 'EEG')
        tUnits = 's';
    elseif strcmp(expmode, 'patch')
        tUnits = 'ms';
    end
end
% TODO: Get rid of chUnits
if nargin < 11 || isempty(chUnits)
    if strcmp(expmode, 'EEG')
        if size(d, 2) == 1
            chUnits = 'mV';
        else
            chUnits = 'uV';
        end
    elseif strcmp(expmode, 'patch')
        [~, chUnits, ~] = identify_channels(d);
    end
elseif ndim == 2 && ~ischar(chUnits) ...
    || ndim == 3 && ~iscellstr(chUnits)
    error('chUnits must be a char array for 2-d data and a cell array of char arrays for 3-d data!');
end
if nargin < 11 || isempty(chLabels)
    if strcmp(expmode, 'EEG')
        chLabels = 'EEG amplitude (mV)';
    elseif strcmp(expmode, 'patch')
        [~, ~, chLabels] = identify_channels(d);
    end
elseif ndim == 2 && ~ischar(chLabels) ...
    || ndim == 3 && ~iscellstr(chLabels)
    error('chLabels must be a char array for 2-d data and a cell array of char arrays for 3-d data!');
end

% Find data parameters
if strcmp(tUnits, 'ms')
    si = sius/1e3;              % sampling interval in ms
elseif strcmp(tUnits, 's')
    si = sius/1e6;              % sampling interval in sec
end
fprintf('Sampling interval = %d %s\n', si, tUnits);
ntps = size(d, 1);              % number of time points (samples)
nchannels = size(d, 2);         % number of channels
fprintf('Number of samples = %d\n', ntps);
fprintf('Number of channels = %d\n', nchannels);
if ndim == 3
    nsweeps = size(d, 3);       % number of sweeps
    fprintf('Number of sweeps = %d\n\n', nsweeps);
end

% Set up time vector
tvec = si*(1:ntps)';        % see units for sampling interval

% Set default xlimits
if nargin < 6 || isempty(left)
    left = 0;
end
if nargin < 7 || isempty(right)
    right = tvec(end);
end
fprintf('Interval to show = [%g %g]\n', left, right);

% Set up labels for each trace
if ndim == 2        % Usually EEG
    traceLabels = cell(1, nchannels);
    for j = 1:nchannels
        traceLabels{j} = ['Channel #', num2str(j)];
    end
elseif ndim == 3    % Usually Patch clamp
    traceLabels = cell(1, nsweeps);
    for k = 1:nsweeps
        traceLabels{k} = ['Sweep #', num2str(k)];
    end
end

% Plot raw data all at once
if ndim == 2 && plotmode == 1        % Could be EEG or patch clamp
    if strcmp(expmode, 'EEG')
        fprintf('Plotting all channels ...\n');
        vvecAll = d;
        figureNum = 1;
        h = plot_traces_abf_bt_helper(figureNum, tvec, vvecAll, left, right, tUnits, nchannels);
        title(sprintf('Data for all channels between %.1f %s and %.1f %s', left, tUnits, right, tUnits));
        ylabel(chLabels);
        legend(traceLabels);
        saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_all.png', left, right)), 'png');    
        hold off;
    elseif strcmp(expmode, 'patch')
        % Plot sweeps
        for j = 1:nchannels
            fprintf('Plotting channel #%d for all sweeps ...\n', j);
            vecAll = squeeze(d(:, j, :));
            figureNum = 1*j;
            h = plot_traces_abf_bt_helper(figureNum, tvec, vecAll, left, right, tUnits, 1);
            title(sprintf('Data for all channels between %.1f %s and %.1f %s', left, tUnits, right, tUnits));
            ylabel(chLabels{j});
            legend(traceLabels);
            saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d_all.png', left, right, j)), 'png');
            hold off;
    %        close(h);
        end
    end
%    close(h);
elseif ndim == 3 && plotmode == 1    % Usually Patch clamp
    % Plot sweeps
    for j = 1:nchannels
        fprintf('Plotting channel #%d for all sweeps ...\n', j);
        vecAll = squeeze(d(:, j, :));
        figureNum = 100*j;
        h = plot_traces_abf_bt_helper(figureNum, tvec, vecAll, left, right, tUnits, nsweeps);
        title(sprintf('Data for all sweeps between %.1f %s and %.1f %s', left, tUnits, right, tUnits));
        ylabel(chLabels{j});
        legend(traceLabels);
        saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d_all.png', left, right, j)), 'png');
        hold off;
%        close(h);
    end
end

% Plot raw data (each channel and/or sweep individually)
if ndim == 2 && plotmode == 2        % Could be EEG or patch clamp
    if strcmp(expmode, 'EEG')
        for j = 1:nchannels
            fprintf('Plotting channel #%d ...\n', j);
            vvec = d(:, j);
            figureNum = 1+j;
            h = plot_traces_abf_bt_helper(figureNum, tvec, vvec, left, right, tUnits, 1);
            title(sprintf('Data for %s between %.1f %s and %.1f %s', traceLabels{j}, left, tUnits, right, tUnits));
            ylabel(chLabels);
            saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d.png', left, right, j)), 'png');
            hold off;
            close(h);
        end
    elseif strcmp(expmode, 'patch')
        for j = 1:nchannels
            fprintf('Plotting channel #%d ...\n', j);
            vvec = d(:, j);
            figureNum = 1+j;
            h = plot_traces_abf_bt_helper(figureNum, tvec, vvec, left, right, tUnits, 1);
            title(sprintf('Data for %s between %.1f %s and %.1f %s', traceLabels{j}, left, tUnits, right, tUnits));
            ylabel(chLabels{j});
            saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d.png', left, right, j)), 'png');
            hold off;
            close(h);
        end
    end
elseif ndim == 3 && plotmode == 2    % Usually Patch clamp        %%% Need to fix this part too
    if strcmp(expmode, 'EEG')
        for j = 1:nchannels
            for k = 1:nsweeps
                fprintf('Plotting channel #%d and sweep #%d ...\n', j, k);
                vec = d(:, j, k);
                figureNum = 100*j+k;
                h = plot_traces_abf_bt_helper(figureNum, tvec, vec, left, right, tUnits, 1);
                title(sprintf('Data for %s between %.1f %s and %.1f %s', ...
                        traceLabels{k}, left, tUnits, right, tUnits));
                ylabel(chLabels);
                saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d_Sweep%d.png', ...
                        left, right, j, k)), 'png');
                hold off;
                close(h);
            end
        end
    elseif strcmp(expmode, 'patch')
        for j = 1:nchannels
            for k = 1:nsweeps
                fprintf('Plotting channel #%d and sweep #%d ...\n', j, k);
                vec = d(:, j, k);
                figureNum = 100*j+k;
                h = plot_traces_abf_bt_helper(figureNum, tvec, vec, left, right, tUnits, 1);
                title(sprintf('Data for %s between %.1f %s and %.1f %s', ...
                    traceLabels{k}, left, tUnits, right, tUnits));
                ylabel(chLabels{j});
                saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d_Sweep%d.png', ...
                    left, right, j, k)), 'png');
                hold off;
                close(h);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h] = plot_traces_abf_bt_helper(fig_num, tvec, vvec, left, right, tUnits, channels_or_sweep_num)
%% Plots traces, sets axes and xlabel
% Usage: [h] = plot_traces_abf_bt_helper(fig_num, tvec, vvec, left, right, tUnits, channels_or_sweep_num)
%    h            - figure
%    fig_num            - figure of number
%    tvec            - time vector for plotting
%    vvec            - compressed data vector
%     left            - time interval start
%    right            - time interval end
%    tUnits            - units of time
%    channels_or_sweep_num    - number of channels if plotmode == 1, 1 if plotmode == 2
% Used by:
%        /media/shareX/brianT/ABF_plotting/plot_traces_abf_bt.m

minimum = min(min(vvec));
maximum = max(max(vvec));
range = maximum - minimum;
h = figure(fig_num);
set(h, 'Visible', 'Off');
clf(h);
for k = 1:channels_or_sweep_num
    plot(tvec, vvec(:, k));    hold on;
end
if range ~= 0
    axis([left, right, minimum - 0.2 * range, maximum + 0.2 * range]);
else
    xlim([left, right]);
end
xlabel(['Time (', tUnits, ')']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
OLD CODE:

%    minimum = min(min(vvecAll));
%    maximum = max(max(vvecAll));
%    range = maximum - minimum;
%    h = figure(1);
%    set(h, 'Visible', 'Off');
%    clf(h);
%    for j = 1:nchannels
%        plot(tvec, vvecAll(:, j));    hold on;
%    end
%    axis([left right minimum-0.2*range maximum+0.2*range]);
%    xlabel(['Time (', tUnits, ')']);
%
%        minimum = min(min(vecAll));
%        maximum = max(max(vecAll));
%        range = maximum - minimum;
%        h = figure(100*j);
%        set(h, 'Visible', 'Off');
%        clf(h);
%        for k = 1:nsweeps
%            plot(tvec, vecAll(:, k));    hold on;
%        end
%        axis([left right minimum-0.2*range maximum+0.2*range]);
%        xlabel(['Time (', tUnits, ')']);
%
%        minimum = min(vvec);
%        maximum = max(vvec);
%        range = maximum - minimum;
%        h = figure(1+j);
%        set(h, 'Visible', 'Off');
%        clf(h);
%        plot(tvec, vvec, 'k');    hold on;
%        axis([left right minimum-0.2*range maximum+0.2*range]);
%        xlabel(['Time (', tUnits, ')']);
%
%            minimum = min(vec);
%            maximum = max(vec);
%            range = maximum - minimum;
%            h = figure(100*j + k);
%            set(h, 'Visible', 'Off');
%            clf(h);
%            plot(tvec, vec, 'k');    hold on;
%            axis([left right minimum-0.2*range maximum+0.2*range]);
%            xlabel(['Time (', tUnits, ')']);
%         % y-axis labels: chLabels{1} is Voltage, chLabels{2} is Current, chLabels{3} is Conductance
%        if j == ind_current    % if the current channel is "Current"
%            ylabel(sprintf('%s (%s)', chLabels{2}, chUnits{2}));        
%        else            %%% what happens if we have conductance channels too?
%            ylabel(sprintf('%s (%s)', chLabels{1}, chUnits{1}));
%        end
%if strcmp(expmode, 'EEG')
%    if ndim == 1
%        chUnits = 'mV';
%    elseif ndim == 2
%        chUnits = {'uV', 'uV'};
%    elseif ndim == 3
%        chUnits = {'uV', 'uV', 'uV'};
%    end
%    chLabels = 'EEG amplitude';
%elseif strcmp(expmode, 'patch')
%    [chUnits, chLabels] = correct_patch_labels(ndim, d);
%end
%


    if ndim == 2        % Usually EEG         %TODO: Not always true, should be fixed
        tUnits = 's';
    elseif ndim == 3    % Usually Patch clamp    %TODO: Not always true, should be fixed
        tUnits = 'ms';
    end
    if ndim == 2        % Usually EEG
        chUnits = 'uV';
    elseif ndim == 3    % Usually Patch clamp
        chUnits = {'mV', 'pA', 'nS'};
    end
    if ndim == 2        % Usually EEG
        chLabels = 'EEG amplitude';
    elseif ndim == 3    % Usually Patch clamp
        chLabels = {'Voltage', 'Current', 'Conductance'};
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

function [chUnits, chLabels] = correct_patch_labels(data)
%% Finds correct units and labels of channels if data is patch clamp
% Usage: [chUnits, chLabels] = correct_patch_labels(data)
%    data        - raw data vector
% Used by:
%        /media/shareX/brianT/ABF_plotting/plot_traces_abf_bt.m

nchannels = size(data, 2);
chUnits = cell(1, nchannels);
chLabels = cell(1, nchannels);
ranges = zeros(1, nchannels);    % stores maximum absolute range for each channel
peak = zeros(1, nchannels);    % maximum for each channel
avgs = zeros(1, nchannels);    % average for each channel
for j = 1:nchannels        % for each channel
    vecAll = squeeze(data(:, j, :));                % all traces in this channel
    ranges(j) = abs(max(max(vecAll)) - min(min(vecAll)));    % maximum range of all traces in this channel
    peak(j) = max(max(vecAll));
end
avgs = mean(mean(data, 1), 3);
if nchannels == 2
    [~, ind_current] = max(ranges);
    ranges(ind_current) = -Inf;
    [~, ind_voltage] = max(ranges);
    chUnits{ind_voltage} = 'mV';
    chUnits{ind_current} = 'pA';
    chLabels{ind_voltage} = 'Voltage';
    chLabels{ind_current} = 'Current';
elseif nchannels == 3
    if abs(range(1) - range(2)) < 5    % if ranges of first two channels are similar enough, assume both are voltage
        chUnits = {'mV', 'mV', 'pA'};         % default: chUnits = {'mV', 'pA', 'nS'};
        chLabels = {'Voltage', 'Voltage', 'Current'};    % default: chLabels = {'Voltage', 'Current', 'Conductance'};
    else                
        [~, ind_conductance] = max(avgs);    % greatest average is presumably "Conductance" (generally positive)
        ranges(ind_conductance) = -Inf;        % all ranges are positive, remove conductance from current detection
        [~, ind_current] = max(ranges);        % the channel with largest range is presumably "Current"
        ranges(ind_current) = -Inf;        % remove current from ranges, leaving only "Voltage"
        [~, ind_voltage] = max(ranges);
        chUnits{ind_voltage} = 'mV';        % Reorder chUnits and chLabels to new indices
        chUnits{ind_current} = 'pA';
        chLabels{ind_voltage} = 'Voltage';
        chLabels{ind_current} = 'Current';
        chLabels{ind_conductance} = 'Conductance';
        if peak(ind_conductance) > 2000    % Conductance usually exceeds 2000 when pico    
            chUnits{ind_conductance} = 'pS';
        else
            chUnits{ind_conductance} = 'nS';
        end
    end
end

        [chUnits, ~] = correct_patch_labels(d);
        [~, chLabels] = correct_patch_labels(d);

ylabel(sprintf('%s (%s)', chLabels, chUnits));
ylabel(sprintf('%s (%s)', chLabels{j}, chUnits{j}));

%}
