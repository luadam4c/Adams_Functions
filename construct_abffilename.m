function abffilename_full = construct_abffilename(filename)
%% Constructs the full file name to an abf file robustly based on filename and display message if doesn't exist
% Usage: abffilename_full = construct_abffilename(filename)
% Arguments:
%	filename	- can be either the full address or must be in current directory
%				.abf is not needed (e.g. 'B20160908_0004')
%
% Used by:
%		cd/plot_traces_abf.m
%		cd/plot_FI.m
% 
% 2017-04-11 Moved from plot_traces_abf.m
% TODO: use construct_fullfilename and rename this construct_and_check_fullfilename; make 'Extension' and optional parameter
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
%%% TODO

%% Create full path to abf file robustly
[filepath, filebase, ~] = fileparts(filename);
abffilename = strcat(filebase, '.abf');
if isempty(filepath)
	filepath = pwd;
end
abffilename_full = fullfile(filepath, abffilename);
fprintf('Full path to abf file: %s\n', abffilename_full);
if exist(abffilename_full, 'file') ~= 2
	error('This abf file doesn''t exist!');
end

