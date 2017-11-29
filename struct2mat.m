function [m, fileName, logHeader, logVariables] = struct2mat (structure, varargin)
%% Saves each variable in a structure as a variable in a MAT-file and create a logHeader and a logVariables
% Usage: [m, fileName, logHeader, logVariables] = struct2mat (structure, varargin)
% Outputs: TODO
% Arguments: TODO
%	structure	- must be a structure
%	fileName	- file name to save as, can include .mat or not	% TODO: make this an optional argument
%	logHeader	- a column cell array of labels for the fields of the structure % TODO: make this an optional argument
%
% Used by:
%
% 2016-11-07 Envisioned
% 2017-05-21 Created
% 

%% Constants
DEFAULT_FILENAME = 'mystruct.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments %% TODO: inputParser
if nargin < 1
	error('Not enough input arguments, type ''help struct2mat'' for usage');
elseif ~isstruct(structure)
	error('First argument must be a structure array!');
end


%% Perform task
% If no file name provided, create a file name
if nargin < 2 || isempty(fileName)
	% Find original structure name, if any
	structName = inputname(1);

	% Create file name
	if ~isempty(structName)		% original structure name provided
		% Use original structure name to for file name
		fileName = [structName, '.mat'];
	else				% no original structure name provided
		% Use default file name
		fileName = DEFAULT_FILENAME;
	end
end

% Replace any original extension in fileName with 'mat'
[filePath, fileBase, fileExt] = fileparts(fileName);
fileName = fullfile(filePath, [fileBase, '.mat']);

% Save fields as individual variables in a MAT-file named by fileName
save(fileName, '-struct', 'structure', '-v7.3');	

% Return MAT-file object
m = matfile(fileName, 'Writable', true);	% MAT-file object for fileName

% Create a log for variables (the field names of the structure)
m.logVariables = fieldnames(structure);

% Create a log for headers 
if nargin < 2 || isempty(logHeader)
	% Use field names as header
	m.logHeader = m.logVariables;
else
	% TODO: check length of logHeader
	m.logHeader = logHeader;
end

