function print_cellstr (cellStr, varargin)
%% Displays the contents of a cell array of character arrays/strings in separate lines
% Usage: print_cellstr (cellStr, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Side Effects:
%       Prints to standard output
% Arguments:    
%       cellStr     - the cell string to print
%                   must be a a cell array of character arrays/strings
%       varargin    - 'FileID': file ID returned by fopen()
%                   must be an integer
%                   default == 1 (standard output)
%
% Requires:
%
% Used by:    
%       cd/print_structure.m
%       cd/m3ha_import_raw_traces.m
%
% File History:
% 2018-06-21 Created by Adam Lu
% TODO: Do not limit cell string to a vector (need nRows, nColumns)
% TODO: Change delimiter
% 

%% Default values for optional arguments
fileIdDefault = 1;                  % standard output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'cellStr', ...                  % the cell string to print
    @(x) assert(iscell(x) && (all(cellfun(@ischar, x)) ...
                || all(cellfun(@isstring, x))), ...
                ['Second input must be a cell array ', ...
                'of strings or character arrays!']));
%    @(x) iscellstr(x) || isstring(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FileId', fileIdDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'integer'}));

% Read from the Input Parser
parse(iP, cellStr, varargin{:});
fileId = iP.Results.FileId;

%% Do the job
% Count the number of strings to print
nStrings = numel(cellStr);

% Do for each string
for iString = 1:nStrings
    fprintf(fileId, '\t''%s''\n', cellStr{iString});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

disp(strjoin(cellStr, '\n'));

%}
