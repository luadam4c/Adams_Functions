function [matches, nMatches] = sscanf_full (str, formatSpec, varargin)
%% Same as sscanf but treats unmatched parts as whitespace (does not stop until end of string)
% Usage: [matches, nMatches] = sscanf_full (str, formatSpec, varargin)
% Example(s):
%       [A, n] = sscanf_full('test(-23.5 -> 230)', '%f')
%       [A, n] = sscanf_full('test(-23.5 -> 230)', '%d')
% Outputs:
%       matches     - matches of pattern from the input text
%                   specified as a column vector or 2-d array
%       nMatches    - Number of actual matches
%                   specified as a positive integer scalar
% Arguments:    
%       str         - input text to scan
%                   must be a string scalar or a character vector
%       formatSpec  - format of input fields
%                   must be a string scalar or a character vector
%       sizeMatches - (opt) size of matches array
%                   must be a numeric 2d array
%
% Requires:
%
% Used by:    
%       /media/adamX/m3ha/data_dclamp/CountSweeps.m
%       /home/Matlab/minEASE/extract_from_minEASE_output_filename.m
%
% File History:
% 2018-07-31 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
sizeMatchesDefault = Inf;

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

% Add required inputs to the Input Parser
addRequired(iP, 'str', ...                  % input text to scan
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(iP, 'formatSpec', ...               % format of input fields
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'sizeMatches', sizeMatchesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, str, formatSpec, varargin{:});
sizeMatches = iP.Results.sizeMatches;

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
