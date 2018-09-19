function abfFullFileName = construct_abffilename(filename)
%% Constructs the full file name to an abf file robustly based on filename and display message if doesn't exist
% Usage: abfFullFileName = construct_abffilename(filename)
% Arguments:
%	filename	- can be either the full address or must be in current directory
%				.abf is not needed (e.g. 'B20160908_0004')
%
% Used by:
%		cd/parse_abf.m
%		cd/plot_FI.m
% 
% 2017-04-11 Moved from plot_traces_abf.m
% TODO: use construct_fullfilename and rename this construct_and_check_fullfilename; 
% TODO: make 'Extension' an optional parameter
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
%%% TODO

%% Create full path to abf file robustly
[filepath, filebase, ~] = fileparts(filename);
abffilename = strcat(filebase, '.abf');
if isempty(filepath)
	filepath = pwd;
end
abfFullFileName = fullfile(filepath, abffilename);
fprintf('Full path to abf file: %s\n', abfFullFileName);
if exist(abfFullFileName, 'file') ~= 2
	fprintf('This abf file doesn''t exist!');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
