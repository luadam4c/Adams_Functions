function plot_FI (filename, alldata, sius, outfolder)
%% From a current injection protocol, detect spikes for each sweep and make an F-I plot
% Usage: plot_FI (filename, alldata, sius, outfolder)
% Arguments:
%       filename    - must be either the full address or must be in current directory
%                   .abf is not needed (e.g. 'B20160908_0004')
%                   Uses abf2load, and if not found, uses abfload
%       alldata     - (opt) full data
%       sius        - (opt) sampling interval in microseconds
%       outfolder   - (opt) the name of the directory that the plots will be placed
%                   must be a character array
%                   default == a subdirectory named by {filename}_traces
%
% Requires:
%       cd/construct_abffilename.m
%       /home/Matlab/Brians_Functions/identify_channels.m
%       /home/Matlab/Brians_Functions/identify_CI.m
%       /home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%
% Used by:
%       cd/plot_all_abfs.m
% 
% 2017-04-04 - Adapted from analyzeCI_20140730
% 2017-04-04 - Used dirr.m to find all abf files
% 2017-04-04 - Updated file names
% 2017-04-06 - BT - Computed and plotted spike frequency
% 2017-04-11 - Moved a copy of identify_channels.m to /home/Matlab/Adams_Functions/
% 2017-04-11 - Now uses filename directly and calls each file from plot_all_abfs_dir.m
% 2017-04-11 - Cleaned code
% 2017-04-11 - Changed the color map to lines
% 2017-04-13 - BT - Marked on plot spike frequency time interval
% 2017-04-13 - Added alldata, sius as optional arguments
% 2017-05-01 - BT - Converted spif_time to use identify_CI.m
% 2017-06-16 - AL - Updated 
% 2018-01-24 - Added isdeployed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
    error('No filename specified!');
elseif ~ischar(filename)
    error('Filename must be a char array in single quotes!');
end

%% Get file parts
[filepath, filebase, ~] = fileparts(filename);

%% Set default outfolder
if nargin < 2
    % Creates full path to abf file robustly
    abffilename_full = construct_abffilename(filename);

    % Load abf file, si is in us
    if exist('abf2load', 'file') == 2
        [alldata, sius] = abf2load(abffilename_full);
    elseif exist('abfload', 'file') == 2
        [alldata, sius] = abfload(abffilename_full);
    end
end
if nargin < 4 || isempty(outfolder) 
    outfolder = fullfile(filepath, strcat(filebase, '_traces'));
end

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
    addpath(fullfile(functionsdirectory, '/Brians_Functions/'));
                                                % for identify_channels.m
end

%% Create outfolder if not already exists
if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
    fprintf('New directory made: %s\n\n', outfolder);
end

% Extract data parameters
nsweeps = size(alldata, 3);     % Number of sweeps
ntimepoints = size(alldata, 1);  % Number of time points

% Find all spikes
is_local_maximum = zeros(size(alldata));
is_local_minimum = zeros(size(alldata));
is_spike = zeros(size(alldata));
for i = 1:nsweeps
    % Finds all local maxima and minima
    for j = 2:ntimepoints-1
        is_local_maximum(j,1,i) = alldata(j,1,i) > alldata(j-1,1,i) && alldata(j,1,i) >= alldata(j+1,1,i);
        is_local_minimum(j,1,i) = alldata(j,1,i) < alldata(j-1,1,i) && alldata(j,1,i) <= alldata(j+1,1,i);
    end

    % Finds all spike peaks 
    % Criteria for a spike: 
    %  (1) Must be a local maximum 10 mV higher than the previous local minimum
    for j = 2:ntimepoints-1
        if is_local_maximum(j,1,i)
            plmin = j-1;    % possible index of previous local minimum
            while ~ (is_local_minimum(plmin,1,i) || plmin == 1) 
                plmin = plmin - 1;
            end
            if plmin > 1
                % Compare rise in membrane potential to thresholds (10 mV)
                is_spike(j,1,i) = alldata(j,1,i) - alldata(plmin,1,i) > 10;
            end
        end
    end
    %  (2) Must be 5 mV higher than the minimum value between the spike and the following spike
    for j = 2:ntimepoints-1
        if is_spike(j,1,i)
            fspike = j+1;     % possible index of following spike
            while ~ (is_spike(fspike,1,i) || fspike == ntimepoints) 
                fspike = fspike + 1;
            end
            is_spike(j,1,i) = alldata(j,1,i) - min(alldata(j:fspike,1,i)) > 5;
        end
    end
end

% Plot all data together
timepoints = ( sius/1000 : sius/1000 : ntimepoints * sius/1000 )'; % vector for timepoints in msec
h = figure(nsweeps + 1);
set(h, 'Visible', 'off');
clf(h);
figname = [filebase, '_all'];
cm = colormap(lines);
for i = 1:nsweeps
    cdata = alldata(:, 1, i);
    plot(timepoints, cdata, 'Color', cm(mod(i, size(cm, 1)) + 1, :), ...
        'Displayname', ['Sweep #', num2str(i)]);
    hold on;
end
title(strrep(figname, '_', '\_'));
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
legend('Location', 'northeast');
saveas(h, fullfile(outfolder, figname), 'png');
close(h);

% Compute spike frequency for each sweep
spike_freq = zeros(1,nsweeps);
iVecs = squeeze(alldata(:, 2, :));
[~, spif_time] = identify_CI(iVecs, sius);    % time range to take spike frequency (10^-4 s)
parfor i = 1:nsweeps
    cdata = alldata(:,1,i);
    spike_indices = find(is_spike(:,1,i));
    spi_within_injection = spike_indices(spike_indices > spif_time(1) & spike_indices < spif_time(2));
    if length(spi_within_injection) > 1
        spike_freq(i) = (length(spi_within_injection) - 1) / ...
                 ((timepoints(spif_time(2)) - timepoints(spif_time(1)))) * 1000;
    else
        spike_freq(i) = 0;
    end
end

% Plot each sweep individually with detected spikes and computed spike frequency    
parfor i = 1:nsweeps
    cdata = alldata(:,1,i);
    spike_indices = find(is_spike(:,1,i));
    figname = [filebase, '_sweep', num2str(i)];

    h = figure(i);
    set(h, 'Visible', 'off');
    clf(h);
    plot(timepoints, cdata, 'k')
    hold on;
    plot(timepoints(spike_indices), cdata(spike_indices), 'xr');
    y = ylim;
    line([timepoints(spif_time(1)) timepoints(spif_time(1))], [y(1) y(2)]);
    line([timepoints(spif_time(2)) timepoints(spif_time(2))], [y(1) y(2)]);
    plot(timepoints(spike_indices(spike_indices > spif_time(1) & spike_indices < spif_time(2))), ...
          cdata(spike_indices(spike_indices > spif_time(1) & spike_indices < spif_time(2))), 'xg');
    if length(spike_indices) > 1
        text(0.6, 0.85, ['Spike Frequency: ' num2str(spike_freq(i)) ' Hz'], 'Units', 'normalized');
    else
        text(0.6, 0.85, ['Spike Frequency: 0 Hz'], 'Units', 'normalized');
    end
    title(strrep(figname, '_', '\_'));
    xlabel('Time (ms)')
    ylabel('Membrane Potential (mV)')
    saveas(h, fullfile(outfolder, figname), 'png');
    close(h);
end

% Plot spike frequency over current injected (F-I plot)
figname_FI = [filebase, '_FI'];
channelTypes = identify_channels(alldata, 'ExpMode', 'patch');
ind_cur = strcmpi('Current', channelTypes);
currents = zeros(1, nsweeps);
parfor i = 1:nsweeps
    currents(i) = mean(alldata(spif_time(1):spif_time(2), ind_cur, i));    %%%TODO: after identify_CI.m is finished, use to find injection range
end
h = figure(999);
set(h, 'Visible', 'off');
clf(h);
plot(currents, spike_freq);
title(['F-I plot for ', strrep(filebase, '_', '\_')]);
xlabel('Current Injected (pA)');
ylabel('Spike Frequency (Hz)');
saveas(h, fullfile(outfolder, figname_FI), 'png');
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
OLD CODE:    TODO: Move any old versions of the code here in case we need it back in the future

[alldata, si] = abf2load(fullfile(infolder, files(f).name)); % si is the sampling interval in usec

timepoints = ( si/1000 : si/1000 : ntimepoints * si/1000 )';

plot_FI(infolder, outfolder)
%        /home/Matlab/Downloaded_Functions/dirr.m
%% Find all .abf files
[~, ~, filenames] = dirr(infolder, '.abf', 'name');
if isempty(filenames)
    fprintf('No abf files in current directory!\n');
    fprintf('Type ''help plot_FI'' for usage\n');
end
nfiles = numel(filenames);
%parfor f = 1:nfiles            % TODO: Switch to parfor when the code is ready
for f = 1:nfiles            % for debug
end

filebase = [strrep(files(f).name, '.abf', ''), '_sweep', num2str(i)];
filebase = strrep(files(f).name, '.abf', '');
title(strrep(filebase, '_', '\_'));
saveas(h, fullfile(outfolder, filebase), 'png');
filebase = [strrep(files(f).name, '.abf', ''), '_spikecurrent'];
title(strrep(filebase, '_', '\_'));
saveas(h, fullfile(outfolder, filebase), 'png');

legend1 = legend('Sweep #1','Sweep #2','Sweep #3','Sweep #4','Sweep #5','Sweep #6','Sweep #7','Sweep #8','Sweep #9','Sweep #10');
set(legend1,'Position',[0.132589285714285 0.419345238095237 0.209821428571428 0.501488095238095]);
%axis([0 10000 -160 40])
%axis([0 4000 -160 40])

            %             flmin = j+1; % possible index of following local minimum
            %             while ~ (is_local_minimum(flmin,1,i) || flmin == ntimepoints) 
            %                 flmin = flmin + 1;
            %             end
            %             flmax = j+1; % possible index of following local maximum
            %             while ~ (is_local_maximum(flmax,1,i) || flmax == ntimepoints) 
            %                 flmax = flmax + 1;
            %             end
            %             if plmin > 1 && flmin < ntimepoints
            %                is_spike(j,1,i) = alldata(j,1,i) - alldata(plmin,1,i) > 20;
            %                is_spike(j,1,i) = alldata(j,1,i) - alldata(plmin,1,i) > 10 && ...
            %                     (alldata(j,1,i) - alldata(flmin,1,i) > 10 || alldata(flmax,1,i) - alldata(j,1,i) < 10);

    hold off;

% For debug
c_is_local_maximum = is_local_maximum(:,:,1);
c_is_local_minimum = is_local_minimum(:,:,1);
c_is_spike = is_spike(:,:,1);
c_data = alldata(:,1,1);

    plot(timepoints, cdata, 'Color', jmap((i * floor(size(jmap,1)/10)), :), ...
        'Displayname', ['Sweep #', num2str(i)]);

jmap = colormap(jet);

ind_cur = find(vcc == 2);

%}
