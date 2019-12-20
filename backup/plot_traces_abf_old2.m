function [d, sius, tvec] = plot_traces_abf (filename, left, right, plotmode, t_units, ch_units, ch_labels)
%% Takes an abf file and plots all traces
% Usage: [d, sius, tvec] = plot_traces_abf (filename, left, right, plotmode, t_units, ch_units, ch_labels)
%	d		- full data
%	sius		- sampling interval in microseconds
%	tvec		- a time vector that can be used to plot things, units are in t_units (see below for default)
%       filename	- must be either the full address or must be in current directory
%				.abf is not needed (e.g. 'B20160908_0004')
%				Uses abf2load, and if not found, uses abfload
%				sampling interval is assumed to be in microseconds
%       left		- (opt) the start of the time interval of interest (in seconds for 2-d data and in ms for 3-d data)
%				default == 0
%       right		- (opt) the end of the time interval of interest (in seconds for 2-d data and in ms for 3-d data)
%				default == tvec(end)
%	plotmode	- (opt) 1: all traces are to be plotted together;
%				2: each trace is to be plotted individually
%				default == 1
%	t_units		- (opt) units for time, must be a character array
%				default == 's' for 2-d data and 'ms' for 3-d data
%	ch_units	- (opt) units for each channel, must be a char array for 2-d data 
%					and a cell array of char arrays for 3-d data
%				default == 'uV' for 2-d data and {'mV', 'pA', 'nS'} for 3-d data)
%	ch_labels	- (opt) labels for each channel, must be a char array for 2-d data
%					and a cell array of char arrays for 3-d data
%				default == 'EEG amplitude' for 2-d data and {'Voltage', 'Current', 'Conductance'} for 3-d data
%
% Requires:
%		/home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
% Used by:
%		/media/shareX/share/Adam/Sample_files_from_Katie/test_sweeps.m
%
% File history
% 20160922 - adapted from plot_traces_abf_EEG
% 20160213 - BT - added labelling detection between current and voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check for required arguments
if nargin < 1
	error('No filename specified!');
elseif ~ischar(filename)
	error('Filename must be a char array in single quotes!');
end

%% Set default plot mode
if nargin < 4 || isempty(plotmode)
	plotmode = 1;
	fprintf('Using plot mode == %d ...', plotmode);
elseif plotmode ~= 1 && plotmode ~= 2
	fprintf('Using plot mode == %d ...', plotmode);
	error('This plot mode is not currently supported!');
end

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
	functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
	functionsdirectory = '/scratch/al4ng/Matlab/';
else
	error('Valid functionsdirectory does not exist!');
end
addpath_custom(fullfile(functionsdirectory, '/Downloaded_Functions/'));	% for abf2load.m

%% Create full path to abf file robustly
[filepath, filebase, ~] = fileparts(filename);
abffilename = strcat(filebase, '.abf');
if isempty(filepath)
	filepath = pwd;
end
filename_full = fullfile(filepath, abffilename);
if exist(filename_full, 'file') ~= 2
	error('This abf file doesn''t exist!');
end
fprintf('Full path to abf file: %s\n', filename_full);

%% Create output folder to put figures
outfolder = fullfile(filepath, strcat(filebase, '_traces'));
if exist(outfolder, 'dir') ~= 7
	mkdir(outfolder);
	fprintf('New directory made: %s\n\n', outfolder);
end

%% Load abf file, si is in us
if exist('abf2load', 'file') == 2
	[d, sius] = abf2load(filename_full);
elseif exist('abfload', 'file') == 2
	[d, sius] = abfload(filename_full);
end

%% Find data dimensions
ndim = ndims(d);	% number of dimensions in data
if ndim > 3
	error('Cannot parse data with more than 3 dimensions!');
end
fprintf('Number of data dimensions =  %d\n\n', ndim);

%% Units and labels
if nargin < 5 || isempty(t_units)
	if ndim == 2		% Usually EEG
		t_units = 's';
	elseif ndim == 3	% Usually Patch clamp
		t_units = 'ms';
	end
elseif ~ischar(t_units)
	error('t_units must be a char array!');
end
if nargin < 6 || isempty(ch_units)
	if ndim == 2		% Usually EEG
		ch_units = 'uV';
	elseif ndim == 3	% Usually Patch clamp
		ch_units = {'mV', 'pA', 'nS'};
	end
elseif ndim == 2 && ~ischar(ch_units) ...
	|| ndim == 3 && ~iscellstr(ch_units)
	error('ch_units must be a char array for 2-d data and a cell array of char arrays for 3-d data!');
end
if nargin < 7 || isempty(ch_labels)
	if ndim == 2		% Usually EEG
		ch_labels = 'EEG amplitude';
	elseif ndim == 3	% Usually Patch clamp
		ch_labels = {'Voltage', 'Current', 'Conductance'};
	end
elseif ndim == 2 && ~ischar(ch_labels) ...
	|| ndim == 3 && ~iscellstr(ch_labels)
	error('ch_labels must be a char array for 2-d data and a cell array of char arrays for 3-d data!');
end

% Find data parameters
if strcmp(t_units, 'ms')
	si = sius/1e3;		% sampling interval in ms
elseif strcmp(t_units, 's')
	si = sius/1e6;		% sampling interval in sec
end
fprintf('Sampling interval = %d %s\n', si, t_units);
ntps = size(d, 1);		% number of time points (samples)
nchannels = size(d, 2);		% number of channels
fprintf('Number of samples = %d\n', ntps);
fprintf('Number of channels = %d\n', nchannels);
if ndim == 3
	nsweeps = size(d, 3);	% number of sweeps
	fprintf('Number of sweeps = %d\n\n', nsweeps);
end

% Set up time vector
tvec = si*(1:ntps)';		% see units for sampling interval

% Set default xlimits
if nargin < 2 || isempty(left)
	left = 0;
end
if nargin < 3 || isempty(right)
	right = tvec(end);
end
fprintf('Interval to show = [%g %g]\n', left, right);

% Set up labels for each trace
if ndim == 2		% Usually EEG
	trace_labels = cell(1, nchannels);
	for j = 1:nchannels
		trace_labels{j} = ['Channel #', num2str(j)];
	end
elseif ndim == 3	% Usually Patch clamp
	trace_labels = cell(1, nsweeps);
	for k = 1:nsweeps
		trace_labels{k} = ['Sweep #', num2str(k)];
	end
end

% Plot raw data all at once
if ndim == 2 && plotmode == 1		% Usually EEG
	fprintf('Plotting all channels ...\n');
	vvec_all = d;
	minimum = min(min(vvec_all));
	maximum = max(max(vvec_all));
	range = maximum - minimum;
	h = figure(1);
%	set(h, 'Visible', 'Off');
	clf(h);
	for j = 1:nchannels
		plot(tvec, vvec_all(:, j));	hold on;
	end
	axis([left right minimum-0.2*range maximum+0.2*range]);
	title(sprintf('Data for all channels between %.1f %s and %.1f %s', left, t_units, right, t_units));
	xlabel(['Time (', t_units, ')']);
	ylabel(sprintf('%s (%s)', ch_labels, ch_units));
	legend(trace_labels);
	saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_all.png', left, right)), 'png');
	hold off;
%	close(h);
elseif ndim == 3 && plotmode == 1	% Usually Patch clamp
	% Find the channel with maximum range and make it "Current"		%%% MAY NEED TO FIX
	ranges = zeros(1, nchannels);	% stores maximum range for each channel
	for j = 1:nchannels		% for each channel
		vec_all = squeeze(d(:, j, :));				% all traces in this channel
		ranges(j) = abs(max(max(vec_all)) - min(min(vec_all)));	% maximum range of all traces in this channel
	end
	[~, ind_current] = max(ranges);	% the channel with largest range is presumably "Current"

	% Plot sweeps
	for j = 1:nchannels
		fprintf('Plotting channel #%d for all sweeps ...\n', j);
		vec_all = squeeze(d(:, j, :));
		minimum = min(min(vec_all));
		maximum = max(max(vec_all));
		range = maximum - minimum;
		h = figure(100*j);
%		set(h, 'Visible', 'Off');
		clf(h);
		for k = 1:nsweeps
			plot(tvec, vec_all(:, k));	hold on;
		end
		axis([left right minimum-0.2*range maximum+0.2*range]);
		title(sprintf('Data for all sweeps between %.1f %s and %.1f %s', left, t_units, right, t_units));
		xlabel(['Time (', t_units, ')']);
		% y-axis labels: ch_labels{1} is Voltage, ch_labels{2} is Current, ch_labels{3} is Conductance
		if j == ind_current	% if the current channel is "Current"
			ylabel(sprintf('%s (%s)', ch_labels{2}, ch_units{2}));		
		else			%%% what happens if we have conductance channels too?
			ylabel(sprintf('%s (%s)', ch_labels{1}, ch_units{1}));
		end
		legend(trace_labels);
		saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d_all.png', left, right, j)), 'png');
		hold off;
%		close(h);
	end
	
end

% Plot raw data (each channel and/or sweep individually)
if ndim == 2 && plotmode == 2		% Usually EEG
	for j = 1:nchannels
		fprintf('Plotting channel #%d ...\n', j);
		vvec = d(:, j);
		minimum = min(vvec);
		maximum = max(vvec);
		range = maximum - minimum;
		h = figure(1+j);
		set(h, 'Visible', 'Off');
		clf(h);
		plot(tvec, vvec, 'k');	hold on;
		axis([left right minimum-0.2*range maximum+0.2*range]);
		title(sprintf('Data for %s between %.1f %s and %.1f %s', trace_labels{j}, left, t_units, right, t_units));
		xlabel(['Time (', t_units, ')']);
		ylabel(sprintf('%s (%s)', ch_labels, ch_units));
		saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d.png', left, right, j)), 'png');
		hold off;
		close(h);
	end
elseif ndim == 3 && plotmode == 2	% Usually Patch clamp		%%% Need to fix this part too
	for j = 1:nchannels
		for k = 1:nsweeps
			fprintf('Plotting channel #%d and sweep #%d ...\n', j, k);
			vec = d(:, j, k);
			minimum = min(vec);
			maximum = max(vec);
			range = maximum - minimum;
			h = figure(100*j + k);
			set(h, 'Visible', 'Off');
			clf(h);
			plot(tvec, vec, 'k');	hold on;
			axis([left right minimum-0.2*range maximum+0.2*range]);
			title(sprintf('Data for %s between %.1f %s and %.1f %s', trace_labels{k}, left, t_units, right, t_units));
			xlabel(['Time (', t_units, ')']);
			ylabel(sprintf('%s (%s)', ch_labels{j}, ch_units{j}));
			saveas(h, fullfile(outfolder, sprintf('%.1f_%.1f_Channel%d_Sweep%d.png', left, right, j, k)), 'png');
			hold off;
			close(h);
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Instead of copying and pasting, should place the parts that exist more than once into functions that get called
function prepare_to_plot()
%%%



