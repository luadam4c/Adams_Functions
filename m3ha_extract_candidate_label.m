function candidateLabels = m3ha_extract_candidate_label (strs, varargin)
%% Extracts the cell name from strings but ignores anything before filesep
% Usage: candidateLabels = m3ha_extract_candidate_label (strs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       m3ha_extract_candidate_label(oscDataPaths)
%       m3ha_extract_candidate_label(oscDataPaths, 'FromBaseName', true)
%       m3ha_extract_candidate_label(oscDataPaths, 'ForceSingleOutput', true)
%
% Outputs:
%       candidateLabels - candidate labels
%                       specified as a character vector
%
% Arguments:
%       strs        - strings
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'FromBaseName': whether to restrict to the base name
%                                       (remove everything before last filesep
%                                        and after last .)
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for extract_substrings()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_substrings.m
%
% Used by:
%       cd/m3ha_plot_figure08.m

% File History:
% 2020-04-19 Modified from m3ha_extract_candidate_label.m
% 

%% Hard-coded parameters
candidateLabelPattern = 'candidateIDs_[0-9,-]*';

%% Default values for optional arguments
fromBaseNameDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'strs', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FromBaseName', fromBaseNameDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, strs, varargin{:});
fromBaseName = iP.Results.FromBaseName;

% Keep unmatched arguments for the extract_substrings() function
otherArguments = iP.Unmatched;

%% Do the job
% Extract the cell names
candidateLabels = extract_substrings(strs, 'FromBaseName', fromBaseName, ...
                                'RegExp', candidateLabelPattern, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
