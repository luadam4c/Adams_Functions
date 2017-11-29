function plot_all_abfs_dir (directory, expmode)
%% Plots all abf files in directory
% Usage: plot_all_abfs_dir (directory, expmode)
%
% Arguments:	directory	- (opt) the name of the directory containing the abf files, e.g. '20161216'
%				must be a character array
%				default == pwd
%		expmode		- (opt)	'EEG'
%					'patch'
%				must be consistent with plot_traces_abf.m
%				default == 'patch'
%
% Requires:
%		cd/plot_traces_abf.m
%		cd/plot_FI.m
%		/home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%		/home/Matlab/Brians_Functions/identify_CI.m
%		
% File history: 
% 2016-09-22 Created
% 2017-04-11 Added expmode as arguments
% 2017-04-11 Now uses dirr.m to find abf files in subdirectories too
% 2017-04-17 - BT - Creates F-I plot for current injection protocols
% 2017-04-19 - BT - Changed detection method to difference of sweep averages

max_swp_spacing = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
	functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
	functionsdirectory = '/scratch/al4ng/Matlab/';
else
	error('Valid functionsdirectory does not exist!');
end
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));	% for dirr.m, abf2load.m or abfload.m
addpath(fullfile(functionsdirectory, '/Brians_Functions/'));		% for identify_CI.m

%% Set defaults
if nargin < 1
	directory = pwd;
end
if nargin < 2
	expmode = 'patch';
end

%% Find all .abf files
[~, ~, filenames] = dirr(directory, '.abf', 'name');
if isempty(filenames)
	fprintf('No abf files in current directory!\n');
	fprintf('Type ''help plot_all_abfs_dir'' for usage\n');
end
nfiles = numel(filenames);

%% Plot traces from each file using plot_traces_abf.m
parfor k = 1:nfiles
	% Plot all traces
	[d, sius] = plot_traces_abf(filenames{k}, expmode);
	[is_CI, CI_range] = identify_CI(d);
	% If it's a current injection protocol, detect spikes for each sweep and make an F-I plot
	if is_CI
		plot_FI(filenames{k}, d, sius);
	end
end


%{
OLD CODE:

files = dir(directory);
if strfind(files.name, '.abf')

files = dirr(directory, '.abf');
for file = files'

	if nargin >= 3
		plot_traces_abf(fullfile(directory, filenames{k}), expmode, recmode);
	else
		plot_traces_abf(fullfile(directory, filenames{k}), expmode);
	end

	% Load abf file
	abffilename_full = construct_abffilename(filenames{k});	% creates full path to abf file robustly
	if exist('abf2load', 'file') == 2
		[d, sius] = abf2load(abffilename_full);
	elseif exist('abfload', 'file') == 2
		[d, sius] = abfload(abffilename_full);
	end

		injection_data = current_data(12000:20000,:,:);
		if std(injection_data,0,1) < 4 & abs(max(max(injection_data)) - min(min(injection_data))) > 100
			plot_FI_bt(filenames{k}, d, sius);
		end
	if length(current_data) > 20000			% tests sweeps if values within injection time range are consistent
		
		injection_data = squeeze(current_data(12000:20000, 1, :));	% current values within typical injection time range
		avgs_byswp = mean(injection_data, 1);				% average of sweeps
		reduction = abs(diff(diff(avgs_byswp)));			% reduces sweeps into differences between successive sweep averages		%%% TODO: this is not differences but rather differences of differences
		[~, max_swp] = max(avgs_byswp);					% highest sweep by average
		max_swp_peaks_avg = mean(findpeaks(injection_data(:, max_swp)));	% average peak value of greatest sweep
		if reduction < max_swp_spacing & max_swp_peaks_avg > 100 & size(d, 3) > 1	% sweep avgs should be separated by constant
			plot_FI(filenames{k}, d, sius);
		end
	end
%}
