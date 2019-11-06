function abfFullFileName = construct_abffilename(filename)
%% Constructs the full file name to an abf file robustly based on filename and display message if doesn't exist
% Usage: abfFullFileName = construct_abffilename(filename)
% Arguments:
%	filename	- can be either the full address or must be in current directory
%				.abf is not needed (e.g. 'B20160908_0004')
%
% Used by:
%		cd/parse_abf.m
%		cd/parse_current_injection_protocol.m

% File History: 
% 2017-04-11 - Moved from plot_traces_abf.m
% TODO: use construct_fullfilename and rename construct_abffilename -> construct_and_check_fullfilename; 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'filename', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B
% Read from the input Parser
parse(iP, filename, varargin{:});

%% Create full path to abf file robustly
% Get the file directory and file base, ignoring the file extension
[fileDir, filebase, ~] = fileparts(filename);

% Append .abf to the file base, then reconstruct the full file name
abffilename = fullfile(fileDir, strcat(filebase, '.abf')];


if isempty(fileDir)
	fileDir = pwd;
end
abfFullFileName = fullfile(fileDir, abffilename);
fprintf('Full path to file: %s\n', abfFullFileName);

% Return an empty string if the file doesn't exist
if exist(abfFullFileName, 'file') ~= 2
    abfFullFileName = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fprintf('Full path to abf file: %s\n', abfFullFileName);
if exist(abfFullFileName, 'file') ~= 2
    fprintf('This abf file doesn''t exist!');
    return;
end

%}
