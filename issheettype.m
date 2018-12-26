function [results, sheettypes] = issheettype (candidates, varargin)
%% Check whether a string or each string in a cell array is a valid spreadsheet type accepted by readtable()
% Usage: [results, sheettypes] = issheettype (candidates, varargin)
% Outputs:    
%       results     - indication of whether the specified string is
%                        a spreadsheet type
%                   specified as a logical array
%       sheettypes  - validated sheettypes, if any
%                   specified as a string vector, a character vector, 
%                       or a cell array of character vectors
%                   returns the shortest match if matchMode == 'substring' 
%                       (sames as validatestring())
% Arguments:
%       candidates  - string or strings to check
%                   must be a string vector, a character vector, 
%                       or a cell array of character vectors
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
%       cd/all_swd_sheets.m
%       cd/atf2sheet.m
%       cd/create_waveform_train.m
%       cd/create_pulse_train_series.m
%       cd/load_params.m
%       cd/load_swd_sheets.m
%       cd/parse_all_abfs.m
%       cd/parse_all_swds.m
%       cd/parse_assyst_swd.m
%       cd/parse_atf_swd.m
%       cd/plot_swd_raster.m
%       cd/save_params.m
%       cd/ZG_extract_all_IEIs.m
%       cd/ZG_extract_IEI_thresholds.m
%       cd/ZG_compute_IEI_thresholds.m
%       /home/Matlab/function_template.m
%       /home/Matlab/EEG_gui/combine_EEG_gui_outputs.m
%
% File History:
% 2018-05-15 Modified from issheettype.m
% 2018-05-16 Now uses istype.m
% 2018-10-21 Now removes any '.' in the string candidates
% 

%% Hard-coded parameters
validSheetTypes = {'txt', 'dat', 'csv', 'xls', 'xlsb', 'xlsm', ...
                    'xlsx', 'xltm', 'xltx', 'ods'};
                                        % accepted by readtable()
                                        % Note: from Matlab 2018a Documentation

%% Default values for optional arguments
validateModeDefault = false;        % don't throw error by default
matchModeDefault = 'substring';     % match as a substring by default

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
addParameter(iP, 'ValidateMode', validateModeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchMode', matchModeDefault, ...
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
[results, sheettypes] = istype(candidates, validSheetTypes, ...
                               'ValidateMode', validateMode, ...
                               'MatchMode', matchMode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Requires:
%       /home/Matlab/Adams_Functions/validate_string.m

%% Check strings and validate
if iscell(strings)
    sheettypes = cellfun(@(x) validate_string(x, validSheetTypes, ...
                                            'ValidateMode', validateMode, ...
                                            'MatchMode', matchMode), ...
                                            strings, 'UniformOutput', false);
    results = ~cellfun(@isempty, sheettypes);
elseif ischar(strings)
    sheettypes = validate_string(strings, validSheetTypes, ...
                                'ValidateMode', validateMode, ...
                                'MatchMode', matchMode);
    results = ~isempty(sheettypes);
else
    error(['input argument is in the wrong format! ', ...
            'Type ''help %s'' for usage'], mfilename);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%