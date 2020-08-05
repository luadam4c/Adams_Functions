function [matches, nMatches] = sscanf_full (str, formatSpec, varargin)
%% Same as sscanf but treats unmatched parts as whitespace (does not stop until end of string)
% Usage: [matches, nMatches] = sscanf_full (str, formatSpec, varargin)
% Example(s):
%       [matches, nMatches] = sscanf_full('test(-23.5 -> 230)', '%f')
%       [matches, nMatches] = sscanf_full('test(-23.5 -> 230)', '%d')
%       [matches, nMatches] = sscanf_full({'sim24', 'sim45'}, '%d')
%       [matches, nMatches] = sscanf_full({'sim24-26', 'sim45'}, '%d')
%
% Outputs:
%       matches     - matches of pattern from the input text
%                   specified as a column vector or 2-d array
%       nMatches    - Number of actual matches
%                   specified as a positive integer scalar
%
% Arguments:    
%       str         - input text to scan
%                   must be a string scalar or a character vector
%       formatSpec  - format of input fields
%                   must be a string scalar or a character vector
%       sizeMatches - (opt) size of matches array
%                   must be a numeric 2d array
%
% Requires:
%       cd/array_fun.m
%
% Used by:
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/parse_atf_swd.m
%       cd/read_data_atf.m
%       /media/adamX/m3ha/data_dclamp/CountSweeps.m
%       /home/Matlab/minEASE/extract_from_minEASE_output_filename.m

% File History:
% 2018-07-31 Created by Adam Lu
% 2020-04-20 Now accepts cell arrays and string arrays as input

%% Hard-coded parameters

%% Default values for optional arguments
sizeMatchesDefault = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'str', ...                  % input text to scan
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'formatSpec', ...               % format of input fields
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'sizeMatches', sizeMatchesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, str, formatSpec, varargin{:});
sizeMatches = iP.Results.sizeMatches;

%% Recursively apply if input is a cell array or a string array
if iscell(str) || isstring(str)
    try
        [matches, nMatches] = ...
            array_fun(@(x) sscanf_full(x, formatSpec, varargin{:}), str);
    catch
        [matches, nMatches] = ...
            array_fun(@(x) sscanf_full(x, formatSpec, varargin{:}), str, ...
                        'UniformOutput', false);
    end
    return
end

%% Preparation
% Get maximum number of matches
maxNMatches = sum(sizeMatches);

% If sizeMatches is not infinite, preallocate
if ~isinf(maxNMatches)
    matches = zeros(maxNMatches, 1);
else
    matches = [];
end

%% Find all matches
% Iteratively call sscanf and remove characters
while length(str) > 0
    % Scan for more matches
    [newMatches, nNewMatches, errorMessage, nextIndex] = ...
        sscanf(str, formatSpec, maxNMatches);

    % If there are new matches, append to previous matches
    matches = [matches; newMatches];

    % Remove the scanned parts of the string
    str = str((nextIndex+1):end);

    % Decrease the maximum number of matches for the next iteration
    maxNMatches = maxNMatches - nNewMatches;
end

% Count the number of matches
nMatches = length(matches);

%% If sizeMatches is a finite 2-element array, 
%   reshape matches array to match it
if length(sizeMatches) > 1 && ~isinf(maxNMatches)
    matches = reshape(matches, sizeMatches);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
