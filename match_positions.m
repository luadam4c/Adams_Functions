function element = match_positions(cellArray, cellStr, strToMatch)
%% Finds element(s) of a cell array that matches the positions of elements of another cell array containing specified string(s)
% Usage: element = match_positions(cellArray, cellStr, strToMatch)
% Explanation:
%       TODO
% Example(s):
%       match_positions({45, 15, 2}, {'cars', 'dogs', 'be'}, {'a', 'o'});
%       label = match_positions(labels, types, type);
% Outputs:
%       element     - matched elements
% Arguments:
%       cellArray   - a cell array
%                   must be a cell array
%       cellStr     - a cell array of strings
%                   must be a cell array of strings/character arrays
%       strToMatch  - string(s) to match in the second cell array
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%
% Used by:
%       cd/extract_channel.m

% File History:
% 2018-12-15 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'cellArray', @iscell);
addRequired(iP, 'cellStr', ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['cellStr must be a cell array of character arrays', ...
                'or a string array!']));
addRequired(iP, 'strToMatch', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['strToMatch must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Read from the Input Parser
parse(iP, cellArray, cellStr, strToMatch);

%% Do the job
% Determine whether each position in cellStr matches one of strToMatch
isMatched = contains(cellStr, strToMatch);

% Get all the strings of cellArray in this position
element = cellArray(isMatched);

% Return as a character array if there is only one element
if iscell(element) && numel(element) == 1
    element = element{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%