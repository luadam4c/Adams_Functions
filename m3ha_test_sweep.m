% m3ha_test_sweep.m
%% Test for info on a sweep of interest
%
% Requires:
%       cd/m3ha_find_lts.m
%       cd/m3ha_plot_traces_mat.m
%
% 2016-09-?? Created
% 2016-10-14 Updated outputs for m3ha_find_lts
% 2016-10-18 Added iterations to test parfor
% 2016-10-19 Added infolder, sweeps
%

%% Input/Output folders
infolder = '/media/adamX/m3ha/data_dclamp/take4/matfiles/';
outfolder = '/media/adamX/m3ha/data_dclamp/take4/test_sweeps/';

%% Sweeps to plot
smallest_action_potential = 'B091810_0006_19';			% 2016-10-19
sp_thr_incorrect = 'D092710_0006_5';				% 2016-10-19
sweep_to_plot = sp_thr_incorrect;

%% Sweeps to test LTS algorithm on
spp_incorrect = {'B091810_0005_9', 'B091810_0006_4', 'B091810_0006_5', ...
		'B091810_0006_9', 'B091810_0006_10', 'B091810_0006_14', ...
		'B091810_0006_19', 'B091810_0006_20', 'B091810_0006_24', ...
		'H101310_0002_6', 'H101310_0002_10'};		% 2016-10-19
sp_thr_incorrect = {'D092710_0006_5'};				% 2016-10-19
sweeps = spp_incorrect;

%% For testing parfor
iterations = 1;
% iterations = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot all traces
filename_full = fullfile(infolder, [sweep_to_plot, '.mat']);	% full path to matfile
[tvec0, gvec0, ivec0, vvec0, ...
	tvec1, gvec1, ivec1, vvec1, ...
	tvec2, gvec2, ivec2, vvec2, ...
	tvec3, gvec3, ivec3, vvec3] = ...
	m3ha_plot_traces_mat (filename_full);

for i = 1:iterations
	tic;
	parfor k = 1:numel(sweeps)
%	for k = 1:numel(sweeps)
		% Import data
		filebase = sweeps{k};					% sweep name to test
		filename_full = fullfile(infolder, [filebase, '.mat']);	% full path to matfile
		m = matfile(filename_full);
		tvec0 = m.d_orig(:, 1);
		tvec2 = m.d_mfrs(:, 1);
		vvec0 = m.d_orig(:, 4);
		vvec1 = m.d_mf(:, 4);
		vvec2 = m.d_mfrs(:, 4);
		vvec3 = m.d_mfmaf(:, 4);

		% Find LTS
		[actVhold, maxnoise, peaktime, peak2ndder, peakprom, peakwidth, ...
			peakclass, spikesperpeak, ltspeaktime, ltspeakval, ...
			maxslopetime, maxslopeval, bursttime, spikesperburst] = ...
			m3ha_find_lts (tvec0, vvec0, 1000, 1100, [200 1000], [1000 7960], 1, ...
					outfolder, filebase, tvec2, vvec1, vvec2, vvec3);
		close all
	end
	toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE

%fn_now = 'A092910_0001_23.mat';
%fn_now = 'F101210_0000_1.mat';
%filename_full = '/media/adamX/m3ha/data_dclamp/take4/matfiles/A100110_0008_18.mat';
%filename_full = '/media/adamX/m3ha/data_dclamp/take4/matfiles/A092910_0001_23.mat';
%filename_full = '/media/adamX/m3ha/data_dclamp/take4/matfiles/A092910_0003_18.mat';

[filepath, filebase, fileextension] = fileparts(filename_full);
matfilename = strcat(filebase, fileextension);
fprintf('Filename of sweep is %s\n', filename_full);

%% Load statistics
load('take4/dclampdatalog_take4.mat')

%% Find index of sweep to test
ind = find_in_strings(matfilename, fnrow);
fprintf('Index of sweep is %d\n', ind);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
