function element = match_positions(array, cellStr, strToMatch, varargin)
%% Finds element(s) of an array that matches the positions of elements of a cell array containing specified string(s)
% Usage: element = match_positions(array, cellStr, strToMatch, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       match_positions({45, 15, 2}, {'cars', 'dogs', 'be'}, {'a', 'o'});
%       label = match_positions(labels, types, type);
%       tau1 = match_positions(coeffValues, coeffNames, 'b');
%       tau2 = match_positions(coeffValues, coeffNames, 'd');
%
% Outputs:
%       element     - matched elements
%
% Arguments:
%       array       - an array
%                   must be an array
%       cellStr     - a cell array of strings
%                   must be a cell array of strings/character arrays
%       strToMatch  - string(s) to match in the cell array
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       varargin    - 'MaxNum': maximum number of positions to match
%                   must be a positive integer scalar or Inf
%                   default == Inf
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_peak_decay.m
%       cd/extract_channel.m

% File History:
% 2018-12-15 Created by Adam Lu
% 2018-12-24 Now accepts any array type as the first element
% 

%% Default values for optional arguments
maxNumDefault = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'array');
addRequired(iP, 'cellStr', ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['cellStr must be a cell array of character arrays ', ...
                'or a string array!']));
addRequired(iP, 'strToMatch', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['strToMatch must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MaxNum', maxNumDefault, ...
    @(x) isinf(x) || ispositiveintegerscalar(x));

% Read from the Input Parser
parse(iP, array, cellStr, strToMatch);
maxNum = iP.Results.MaxNum;

%% Do the job
% Determine whether each position in cellStr matches one of strToMatch
isMatched = contains(cellStr, strToMatch);

% Get all the elements of array in this position
element = array(isMatched);

% Restrict to the maximum number of elements
if numel(array) > maxNum
    element = array(1:maxNum);
end

% If a cell array, return as an element if there is only one element
if iscell(element) && numel(element) == 1
    element = element{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%