function m3ha_append_lts_properties (fullmatfilepath)
%% Generates vectors of peak features restricted to those with LTS
% Usage: m3ha_append_lts_properties (fullmatfilepath)
% Arguments: 
%		fullmatfilepath	- (opt) full path to matfile to alter
%				default: //media/adamX/m3ha/data_dclamp/take4/dclampdatalog_take4.mat
%		
% Used by:	
%		cd/m3ha_parse_dclamp_data.m

% File History:
% 2016-10-31 Created
% 2017-12-27 Fixed logheader and logvariables

%% Specify which sweep info to restrict
oldvariables = {'peaktime', 'peak2ndder', 'peakprom', 'peakwidth'};

%% New headers and variables
newheader = {'LTS peak time (ms)', 'LTS peak 2nd derivative (V^2/s^2)', ...
	        'LTS peak prominence (mV)', 'LTS peak width at half-prom (ms)'};
newvariables = {'ltspeaktime', 'ltspeak2ndder', 'ltspeakprom', 'ltspeakwidth'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin >= 1 && exist(fullmatfilepath, 'file') ~= 2
	error('matfile does not exist!');
end

%% Set defaults for optional arguments
if nargin < 1
	fullmatfilepath = '/media/adamX/m3ha/data_dclamp/take4/dclampdatalog_take4.mat';
end

%% Import matfile and make writable
m = matfile(fullmatfilepath, 'Writable', true);
ltspeaktime = m.ltspeaktime;
for v = 1:numel(oldvariables)
	% Load old variable
	temp = m.(oldvariables{v});

	% Modify variable
	numswps = numel(temp);
	for k = 1:numswps
		if isnan(ltspeaktime(k))
			temp(k) = NaN;
		end
	end

	% Check ltspeaktime for sanity
	if strcmp(newvariables{v}, 'ltspeaktime')
		for k = 1:numswps
			if isnan(ltspeaktime(k)) && ~isnan(temp(k)) ...
				|| ~isnan(ltspeaktime(k)) && isnan(temp(k)) ...
				|| ~isnan(ltspeaktime(k)) && ~isnan(temp(k)) && ltspeaktime(k) ~= temp(k)
				fprintf('Old: %g ~= New: %g\n', ltspeaktime(k), temp(k));
				error('Something is wrong in the dclampdatalog file!');
			end
		end
		m.ltsonsettime = temp;		% Save "LTS onset time" as a separate instance
	end	
	
	% Save new variable
	m.(newvariables{v}) = temp;
end

%% Update logheader and logvariables
m.logheader = [m.logheader, newheader];
m.logvariables = [m.logvariables, newvariables];

%{
%% OLD CODE

% This will create problems if m3ha_append_lts_properties is ran more than once
m.logheader = {m.logheader, newheaderentries};
m.logheader = {m.logvariables, newvariables};

%%% DOESN'T WORK
%% Set matfile to be nonwritable
set(m, 'Writable', false);

newheaderentries = {'LTS peak time (ms)', 'LTS peak 2nd derivative (V^2/s^2)', ...
		'LTS peak prominence (mV)', 'LTS peak width at half-prom (ms)'};

%}

