function varargout = find_in_strings (cand, strList, varargin)
%% Returns all indices of a particular string (could be represented by substrings) in a list of strings
% Usage: [indices, matched] = find_in_strings (cand, strList, varargin)
% Explanation:
%   There are three main search modes (parameter 'SearchMode'):
%       'substrings': allows the candidate to be a substring or substrings 
%                       of a match in strList.
%       'exact': candidate must be an exact match
%       'regexp': candidate is a regular expression
%   The latter two cases are similar to strcmp()/strcmpi() or regexp()/regexpi()
%   However, find_in_strings returns indices instead of logical arrays,
%       and optionally returns the matched elements as the second output.
%   Use ismember_custom.m to simply test whether cand is in strList
%   Use find_first_match.m to treat each element of cand
%       as a different candidate
%
% Example(s):
%       strs1 = {'Mark''s fish', 'Peter''s fish', 'Katie''s sealion'};
%       strs2 = ["Mark's fish", "Peter's fish", "Katie's sealion"];
%       find_in_strings('fish', strs1)
%       find_in_strings('Peter', strs2)
%       find_in_strings({'Katie', 'lion'}, strs2)
%       find_in_strings("fish", strs1, 'MaxNum', 1)
%       find_in_strings("Fish", strs1, 'IgnoreCase', 1)
%       find_in_strings('Fish', strs2, 'IgnoreCase', false)
%       find_in_strings("sealion", strs1, 'SearchMode', 'exact')
%       find_in_strings('sea', strs2, 'SearchMode', 'exact', 'ReturnNaN', true)
%       find_in_strings("sea.*", strs1, 'SearchMode', 'reg')
%       find_in_strings('sea.*', strs2, 'SearchMode', 'reg')
%
% Outputs:
%       indices     - indices of strList matching the candidate
%                       could be empty (or NaN if 'ReturnNan' is true)
%                   specified as a positive integer array
%       matched     - matched elements of strList
%                   specified as a cell array if more than one indices 
%                       or the element if only one index; or an empty string
% Arguments:
%       cand        - candidate string or substring(s)
%                       If cand is a list of substrings, all substrings must 
%                           exist in the string to be matched
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       strList     - a list of strings
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       varargin    - 'SearchMode': the search mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'         - cand must be identical to 
%                                           an element in strList
%                       'substrings'    - cand can be a substring or 
%                                           a list of substrings
%                       'regexp'        - cand is considered a regular expression
%                   if search mode is 'exact' or 'regexp', 
%                       cand cannot have more than one elements
%                   default == 'substrings'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxNum': maximum number of indices to find
%                   must be empty or a positive integer scalar
%                   default == numel(strList)
%                   - 'ReturnNan': Return NaN instead of empty if nothing found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/ismatch.m
%
% Used by:
%       cd/char2rgb.m
%       cd/create_pleth_EEG_movie.m
%       cd/ispositiveintegerscalar.m
%       cd/increment_editbox.m
%       cd/m3ha_correct_unbalanced_bridge.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_iox.m
%       cd/parse_spike2_mat.m
%       cd/plot_swd_raster.m
%       cd/renamevars.m
%       cd/update_params.m
%       cd/validate_string.m
%       cd/xolotl_compartment_index.m
%       cd/ZG_extract_all_IEIs.m
%       cd/ZG_extract_all_data.m
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/minEASE/combine_eventInfo.m
%       /home/Matlab/minEASE/extract_from_minEASE_output_filename.m
%       /home/Matlab/minEASE/read_params.m
%       /home/Matlab/minEASE/gui_examine_events.m
%       /home/Matlab/EEG_gui/EEG_gui.m
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m
%       /media/adamX/m3ha/data_dclamp/dclampPassiveFitter.m
%       /media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%       /media/adamX/m3ha/data_dclamp/test_sweep.m
%       /media/adamX/m3ha/data_dclamp/remove_E092810_0000.m
%       /media/adamX/m3ha/data_dclamp/compare_statistics.m
%       /media/adamX/m3ha/optimizer4gabab/optimizergui_4compgabab.m
%       /media/adamX/m3ha/optimizer4gabab/optimizer_4compgabab.m
%       /media/adamX/m3ha/optimizer4gabab/compare_neuronparams.m
%       /media/adamX/RTCl/neuronlaunch.m
%       /media/adamX/RTCl/raster_plot.m
%       /media/adamX/RTCl/tuning_curves.m
%       /media/adamX/RTCl/single_neuron.m

% File History:
% 2016-09--- Created
% 2016-10-13 moved to Adams_Functions
% 2016-11-30 Added searchMode
% 2017-04-05 Fixed the size of str_cell so that it can take column or row arrays
% 2017-04-26 Now cand can be a cell array of substrings too
% 2017-04-27 Improved inputParser scheme
% 2017-05-09 Added elements as output
% 2017-05-25 Changed line width and indentation
% 2017-06-09 Fixed the returned element to be of original case
% 2018-05-01 Added MaxNum as a parameter
% 2018-08-02 Added 'regexp' as a SearchMode
% 2019-01-04 Now uses compute_combined_trace.m instead of intersect_over_cells.m
% 2019-01-04 Simplified code with contains()
% 2019-01-09 Added 'ReturnNan' as an optional argument
% 2019-01-09 Now uses find_custom.m
% 2019-01-10 Renamed find_ind_str_in_cell -> find_in_strings
% 2019-01-10 Updated Explanation section
% 2019-01-10 Now fully supports strings in double quotes
% 2019-01-10 Moved code to is_matching_string.m

%% Hard-coded constants
validSearchModes = {'exact', 'substrings', 'regexp'};

%% Default values for optional arguments
searchModeDefault = 'substrings';       % default search mode
ignoreCaseDefault = false;              % whether to ignore case by default
maxNumDefault = [];                     % will be changed to numel(strList)
returnNanDefault = false;   % whether to return NaN instead of empty 
                            %   if nothing found by default

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
addRequired(iP, 'cand', ...             % a string/substrings of interest
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['cand must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'strList', ...          % a list of strings
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strList must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validSearchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNum', maxNumDefault, ...       % maximum number of indices
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'MaxNum must be either empty or a positive integer scalar!'));
addParameter(iP, 'ReturnNan', returnNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, cand, strList, varargin{:});
searchMode = validatestring(iP.Results.SearchMode, validSearchModes);
ignoreCase = iP.Results.IgnoreCase;
maxNum = iP.Results.MaxNum;
returnNan = iP.Results.ReturnNan;

%% Preparation
% Check relationships between arguments
if ~ischar(cand) && numel(cand) > 1 && ...
    (strcmp(searchMode, 'exact') || strcmp(searchMode, 'regexp'))
    error(['First input cannot have more than one members if ', ...
            '''SearchMode'' is ''exact'' or ''regexp''!']);
end

% Translate search mode to match mode
if strcmp(searchMode, 'substrings')
    matchMode = 'parts';
else
    matchMode = searchMode;
end

%% Use ismatch.m
if nargout >= 2
    [~, indices, matched] = ismatch(strList, cand, ...
                        'MatchMode', matchMode, 'IgnoreCase', ignoreCase, ...
                        'MaxNum', maxNum, 'ReturnNan', returnNan);
else
    [~, indices] = ismatch(strList, cand, ...
                        'MatchMode', matchMode, 'IgnoreCase', ignoreCase, ...
                        'MaxNum', maxNum, 'ReturnNan', returnNan);
end

%% Outputs
varargout{1} = indices;
if nargout >= 2
    varargout{2} = matched;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
