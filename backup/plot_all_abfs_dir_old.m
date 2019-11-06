function plot_all_abfs_dir (directory, expmode)
%% Plots all abf files in directory
% Usage: plot_all_abfs_dir (directory, expmode)
%
% Arguments:
%       directory	- (opt) the name of the directory containing the abf files, e.g. '20161216'
%			must be a character array
%			default == pwd
%	expmode		- (opt)	'EEG'
%				'patch'
%			must be consistent with plot_traces_abf.m
%			default == 'patch'
%
% Requires:
%		cd/plot_traces_abf.m
%		/home/Matlab/Downloaded_Functions/dirr.m
%
% File history: 
% 2016-09-22 Created
% 2017-04-11 Added expmode as arguments
% 2017-04-11 Now uses dirr.m to find abf files in subdirectories too
%

% Set defaults
if nargin < 1
	directory = pwd;
end
if nargin < 2
	expmode = 'patch';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
	functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
	functionsdirectory = '/scratch/al4ng/Matlab/';
else
	error('Valid functionsdirectory does not exist!');
end
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));	% for dirr.m

%% Find all .abf files
[~, ~, filenames] = dirr(directory, '.abf', 'name');
if isempty(filenames)
	fprintf('No abf files in current directory!\n');
	fprintf('Type ''help plot_all_abfs_dir'' for usage\n');
end
nfiles = numel(filenames);

%% Plot traces from each file using plot_traces_abf.m
parfor k = 1:nfiles
	plot_traces_abf(filenames{k}, expmode);
	% TODO: If it's a current injection protocol, detect spikes for each sweep and make an F-I plot
	% if XXX
	%	parse_current_family(filenames{k});
	% end
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

%}
