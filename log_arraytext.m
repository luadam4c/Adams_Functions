function fid = log_arraytext (filename, array, varargin)
%% Create a text file that logs the array information
% Usage: fid = log_arraytext (filename, array, varargin)
% Outputs:
%       fid         - file ID returned by fopen; -1 if file cannot be opened
%                   must be an integer scalar
% Arguments    
%       filename    - file name of the text file for logging
%                   must be a string scalar or a character vector
%       array       - array to be logged
%                   must be a numeric vector, a cell vector, 
%                       a string vector or a struct vector
%
% Requires:
%       cd/check_dir.m
%
% Used by:    
%        ~/m3ha/optimizer4gabab/optimizergui_4compgabab.m
%
% File History:
% 2017-05-01 Created by Adam Lu
% 2018-11-15 Cleaned code; now uses check_dir.m
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'filename', ... % file name of the text file for logging
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'array', ...    % array to be logged
    @(x) validateattributes(x, {'numeric', 'cell', 'string', 'struct'}, ...
                               {'vector'}));

% Read from the input Parser
parse(iP, filename, array, varargin{:});

%% Preparation
% Get the output directory name from the file name
[outfolder, ~, ~] = fileparts(filename);

% Check if needed output directory exist, otherwise create it
check_dir(outfolder);

%% Do the job
% Open file for writing
fid = fopen(filename, 'w');

% Log array according to type
if iscell(array)
    % Display each element in turn
    for k = 1:numel(array)
        if ischar(array{k}) || isstring(array{k})
            % Print the string on its own line
            fprintf(fid, '%s\n', array{k});
        elseif isnumeric(array{k}) && isscalar(array{k})
            % Print the number on its own line
            fprintf(fid, '%g\n', array{k});
        else
            % Display the element
            disp(array{k});
        end
    end
else
    % Display the array
    disp(array)
end

% Close file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
    fprintf('New directory made: %s\n\n', outfolder);
end

if isnumeric(array)
% 'numeric' vector
fprintf(fid, '%g\n', array(k));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%