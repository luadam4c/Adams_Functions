function [fid] = log_arraytext (filename, array, varargin)
%% Create a text file that logs the array information
% Usage: [fid] = log_arraytext (filename, array, varargin)
% Outputs:	fid		- file ID returned by fopen; -1 if file cannot be opened
% Arguments:	
%       filename	- file name of the text file for logging
%				    must be a string scalar or a character vector
%		array		- array to be logged
%				    must be a numeric vector, a cell vector, 
%                       a string vector or a struct vector
%
% Used by:	
%		/media/adamX/m3ha/optimizer4gabab/optimizergui_4compgabab.m
%
% File History:
% 2017-05-01 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
	error('Not enough input arguments, type ''help log_arraytext'' for usage');
end

% Add required inputs to an input Parser
iP = inputParser;
addRequired(iP, 'filename', ... % file name of the text file for logging
	@(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%	@(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'array', ...    % array to be logged
	@(x) validateattributes(x, {'numeric', 'cell', 'string', 'struct'}, ...
	                           {'vector'}));

% Read from the input Parser
parse(iP, filename, array, varargin{:});

%% Check if needed output directory exist, otherwise create it
[outfolder, ~, ~] = fileparts(filename);
if exist(outfolder, 'dir') ~= 7
	mkdir(outfolder);
	fprintf('New directory made: %s\n\n', outfolder);
end

%% Open file to log array according to type
fid = fopen(filename, 'w');
if isnumeric(array)		% 'numeric' vector
	fprintf(fid, '%g\n', array(k));
else				% 'cell', 'string' or 'struct' vector
	for k = 1:numel(array)
		if ischar(array{k}) || isstring(array{k})
			fprintf(fid, '%s\n', array{k});
		elseif isnumeric(array{k}) && isscalar(array{k})
			fprintf(fid, '%g\n', array{k});
		else
			disp(array{k});
		end
	end
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
