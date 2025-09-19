function [results, linestyles] = islinestyle (candidates, varargin)
%% Check whether a string or each string in a cell array is a valid line style accepted by plot() or line()
% Usage: [results, linestyles] = islinestyle (candidates, varargin)
% Outputs:    
%       results     - indication of whether the specified string is
%                        valid linestyle accepted by plot() or line()
%                   specified as a logical array
%       linestyles  - validated linestyles, if any
%                   specified as a string vector, a character vector, 
%                       or a cell array of character vectors
%                   returns the shortest match if matchMode == 'substring' 
%                       (sames as validatestring())
% Arguments:
%       candidates  - string or strings to check
%                   must be a string vector, a character vector, 
%                       or a cell array of character vectors
%       varargin    - 'ValidateMode': whether to validate string and 
%                       throw error if string is not a substring of a line style
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match 
%                       to one of the following:
%                       'exact'         - string must be exact
%                       'substring'     - string can be a substring
%                   if 'ValidateMode' is 'true', matching mode is 
%                       automatically 'substring'
%                   default == 'substring'
%
% Requires:
%       cd/istype.m
%
% Used by:
%       /home/Matlab/EEG_gui/combine_EEG_gui_outputs.m
%       cd/plot_regression_line.m
%       cd/plot_traces.m
%       cd/plot_vertical_shade.m
%       cd/plot_window_boundaries.m

% File History:
% 2018-05-16 Modified from issheettype.m
% 2019-08-27 Added empty option
% 

%% Hard-coded parameters
validLineStyles = {'', '-', '--', ':', '-.', 'none'};
                                        % accepted by line() or plot()
                                        % Note: from Matlab 2018a Documentation

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

% Add required inputs to an Input Parser
addRequired(iP, 'candidates', ...               % string or strings to check
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['candidates must be either a string array, ', ...
                'a character array or a cell array of character vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValidateMode', false, ...     % whether to validate string
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchMode', 'substring', ...  % the matching mode
    @(x) any(validatestring(x, {'exact', 'substring'})));

% Read from the Input Parser
parse(iP, candidates, varargin{:});
validateMode = iP.Results.ValidateMode;
matchMode = iP.Results.MatchMode;

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

%% Check candidates and validate with istype.m
[results, linestyles] = istype(candidates, validLineStyles, ...
                               'ValidateMode', validateMode, ...
                               'MatchMode', matchMode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%