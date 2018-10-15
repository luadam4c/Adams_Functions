function [results, legendLocations] = islegendlocation (strings, varargin)
%% Check whether a string or each string in a cell array is a valid legend location or 'suppress' or 'auto'
% Usage: [results, legendLocations] = islegendlocation (strings, varargin)
% Outputs:    
%       results     - indication of whether the specified string is
%                        valid legend location accepted by legend()
%                   specified as a logical array
%       legendLocations  - validated legendLocations, if any
%                   specified as a string/char-vec or 
%                       a cell array of strings/char-vecs
%                   returns the shortest match if matchMode == 'substring' 
%                       (sames as validatestring())
% Arguments:
%       strings     - string or strings to check
%                   must be a string/char-vec or 
%                       a cell array of strings/char-vecs
%       varargin    - 'ValidateMode': whether to validate string and 
%                       throw error if string is not a substring of a sheettype
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
%       cd/plot_cfit_pulse_response.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m

% File History:
% 2018-10-12 Modified from islinestyle.m
% 

%% Hard-coded parameters
validLegendLocations = {'auto', 'suppress', ...
                        'north', 'south', 'east', 'west', ...
                        'northeast', 'northwest', 'southeast', 'southwest', ...
                        'northoutside', 'southoutside', 'eastoutside', ...
                        'westoutside', 'northeastoutside', ...
                        'northwestoutside', 'southeastoutside', ...
                        'southwestoutside', 'best', 'bestoutside', 'none'};
                                        % accepted by legend()
                                        % Note: from Matlab 2018a Documentation

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to an Input Parser
addRequired(iP, 'strings', ...                  % string or strings to check
    @(x) assert(ischar(x) || ...
                iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x) , ...
                ['strings must be either a string/character array ', ...
                'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValidateMode', false, ...     % whether to validate string
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchMode', 'substring', ...  % the matching mode
    @(x) any(validatestring(x, {'exact', 'substring'})));

% Read from the Input Parser
parse(iP, strings, varargin{:});
validateMode = iP.Results.ValidateMode;
matchMode = iP.Results.MatchMode;

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

%% Check strings and validate with istype.m
[results, legendLocations] = istype(strings, validLegendLocations, ...
                               'ValidateMode', validateMode, ...
                               'MatchMode', matchMode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

validLegendLocations = {'', 'suppress', ...
                        'north', 'south', 'east', 'west', ...
                        'northeast', 'northwest', 'southeast', 'southwest', ...
                        'northoutside', 'southoutside', 'eastoutside', ...
                        'westoutside', 'northeastoutside', ...
                        'northwestoutside', 'southeastoutside', ...
                        'southwestoutside', 'best', 'bestoutside', 'none'};
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%