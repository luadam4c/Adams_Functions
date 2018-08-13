function [fields, ntps, nsweeps] = plotEEGmatfile (matfilename, left, right)
%% Plots EEG data from a matfile
% Usage: [fields, ntps, nsweeps] = plotEEGmatfile (matfilename, left, right)
%       filename is a matfile in current directory
% 		Must have the following structure:
% 		The first vector is the time vector t in ms
% 		The rest are voltage vectors for the different channels
%       left is the start of the time interval of interest in seconds
%       right is the end of the time interval of interest in seconds

% File history
% 2016-09-19 - created
% 2017-05-21 - I think this is outdated, used plot_trace_mat.m instead

%% Check for arguments
if nargin < 1
	error('No filename!')
end

%% Create output directory
outfolder = strrep(matfilename, '.mat', '');
if exist(outfolder) ~= 7
	mkdir(outfolder);
end

%% Construct full file name
if isempty(strfind(matfilename, '/'))
	current_folder = pwd;
	matfilename_full = fullfile(current_folder, matfilename);
else
	matfilename_full = matfilename;
end

%% Load mat file
% Must have the following structure:
% The first vector is the time vector t in ms
% The rest are voltage vectors for the different channels
matdata = load(matfilename_full);

%% Find data parameters
fields = fieldnames(matdata);
ntps = length(matdata.t);		% Number of time points
nvectors = numel(fields);		% Total number of vectors
nsweeps = nvectors - 1;			% Total number of channels

%% Set up vector for timepoints in seconds
tps = matdata.t * 1e-3;
if nargin < 2
	left = tps(1);
end
if nargin < 3
	right = tps(end);
end

%% Plot raw data (each sweep individually)
for i = 1:nsweeps
	cdata = matdata.(fields{i+1});
	figure(i)
	plot(tps, cdata, 'k')
	hold on
	axis([left right min(cdata)*1.2 max(cdata)*1.2])
	title(sprintf('Data between %.1f s and %.1f s', left, right));
	xlabel('Time (s)')
	ylabel('EEG amplitude (uV)')
	if nsweeps > 1
		saveas(gcf, fullfile(outfolder, sprintf('/%.1f_%.1f_Sweep#%d.png', left, right, i)), 'png')
	else
		saveas(gcf, fullfile(outfolder, sprintf('/%.1f_%.1f.png', left, right, i)), 'png')
	end
	hold off
end
