function [results, figtypes] = isfigtype (candidates, varargin)
%% Check whether a string or each string in a cell array is a valid figure type accepted by saveas()
% Usage: [results, figtypes] = isfigtype (candidates, varargin)
% Outputs:    
%       results     - indication of whether the specified string is a figure type
%                   specified as a logical array
%       figtypes    - validated figtypes, if any
%                   specified as a string vector, a character vector, 
%                       or a cell array of character vectors
%                   returns the shortest match if matchMode == 'substring' 
%                       (sames as validatestring())
% Arguments:        
%       candidates  - string or strings to check
%                   must be a string vector, a character vector, 
%                       or a cell array of character vectors
%       varargin    - 'ValidateMode': whether to validate string and 
%                       throw error if string is not a substring of a figtype
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
%       cd/compute_and_plot_average_response.m
%       cd/compute_and_plot_all_responses.m
%       cd/create_waveform_train.m
%       cd/create_pulse_train_series.m
%       cd/plot_all_abfs.m
%       cd/plot_bar.m
%       cd/plot_struct.m
%       cd/plot_traces.m
%       cd/plot_traces_abf.m
%       cd/plot_tuning_curve.m
%       cd/plot_tuning_map.m
%       /home/Matlab/minEASE/compute_plot_average_PSC_traces.m
%       /home/Matlab/minEASE/detect_gapfree_events.m
%       /media/adamX/RTCl/raster_plot.m
%       /media/adamX/RTCl/single_neuron.m
%       /media/adamX/RTCl/tuning_curves.m
%       /media/adamX/RTCl/tuning_maps.m
%       /media/adamX/Settings_Matlab/function_template.m
%
% File History:
% 2017-05-09 Created by Adam Lu
% 2017-05-23 Modified linewidth and indentation
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-05-15 possible_figtypes -> validFigTypes
% 2018-05-16 Now uses istype.m
% 2018-10-21 Now removes any '.' in the string candidates
% 

%% Hard-coded parameters
validFigTypes = {'png', 'fig', 'm', 'mfig', ...
                'jpeg', 'tiff', 'tiffn', 'meta', ...
                'bmpmono', 'bmp', 'bmp16m', 'bmp256', ...
                'hdf', 'pbm', 'pbmraw', 'pcxmono', ...
                'pcx24b', 'pcx256', 'pcx16', 'pgm', ...
                'pgmraw', 'ppm', 'ppmraw'};
                                        % accepted by saveas()
                                        % Note: from Matlab 2017a Documentation

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

%% Preparation
% Remove any '.'
candidates = replace(candidates, '.', '');

%% Check candidates and validate with istype.m
[results, figtypes] = istype(candidates, validFigTypes, ...
                             'ValidateMode', validateMode, ...
                             'MatchMode', matchMode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
