function [indices, elements] = find_in_strings(str, cellarray, varargin)
%% Find all indices of a particular string in a cell array
% Usage: [indices, elements] = find_in_strings(str, cellarray, varargin)
% Outputs:
%       indices     - indices of the cell array containing that exact string
%                       or containing a substring or all substrings provided; 
%                       could be empty
%                   specified as a numeric array
%       elements    - elements of the cell array corresponding to those indices
%                   specified as a cell array if more than one indices 
%                       or the element if only one index; or an empty string
% Arguments:
%       str         - a string or a substring or 
%                       a cell array of substrings of interest
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       cellarray   - a cell array that contains strings
%                   must be a cell array of strings/character arrays
%       varargin    - 'SearchMode': the search mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'         - str must be identical to 
%                                           an element in cellarray
%                       'substrings'    - str can be a substring or 
%                                           a cell array of substrings
%                   if searchmode is 'exact', str cannot be a cell array
%                   default == 'substrings'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       /home/Matlab/Adams_Functions/intersectm.m
% Used by:    
%       /media/adamX/m3ha/data_dclamp/dclampPassiveFitter.m
%       /media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%       /media/adamX/m3ha/data_dclamp/test_sweep.m
%       /media/adamX/m3ha/data_dclamp/remove_E092810_0000.m
%       /media/adamX/m3ha/data_dclamp/compare_statistics.m
%       cd/m3ha_optimizergui_4compgabab.m
%       cd/m3ha_optimizer_4compgabab.m
%       /media/adamX/RTCl/neuronlaunch.m
%       /media/adamX/RTCl/m3ha_network_raster_plot.m
%       /media/adamX/RTCl/tuning_curves.m
%       /media/adamX/RTCl/single_neuron.m
%       /home/Matlab/Adams_Functions/update_params.m
%       /home/Matlab/Adams_Functions/isfigtype.m
%       /home/Matlab/Brians_Functions/minEASE/ExceltoGUIConverter.m
%
% 2016-09--- Created
% 2016-10-13 moved to Adams_Functions
% 2016-11-30 Added searchmode
% 2017-04-05 Fixed the size of str_cell so that it can take column or row arrays
% 2017-04-26 Now str can be a cell array of substrings too
% 2017-04-27 Improved inputParser scheme
% 2017-05-09 Added elements as output
% 2017-05-25 Changed line width and indentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check number of arguments (better error message than inputParser)
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help find_in_strings'' for usage']);
end

%% Add required inputs to an input Parser
iP = inputParser;
addRequired(iP, 'str', ...              % a string/substrings of interest
    @(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x), ...
                ['First input must be either a string/character array ', ...
                'or a cell array of strings/character arrays!']));
addRequired(iP, 'cellarray', ...        % a cell array that contains strings
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))), ...
                ['Second input must be a cell array ', ...
                'of strings/character arrays!']));

%% Add parameter-value pairs to the input Parser
addParameter(iP, 'SearchMode', 'substrings', ...    % the search mode
    @(x) any(validatestring(x, {'exact', 'substrings'})));
addParameter(iP, 'IgnoreCase', false, ...           % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

%% Read parameter values from the input Parser
parse(iP, str, cellarray, varargin{:});
searchmode = validatestring(iP.Results.SearchMode, {'exact', 'substrings'});
ignorecase = iP.Results.IgnoreCase;

%% Check more arguments
if strcmpi(searchmode, 'exact') && iscell(str)
    error('First input cannot be a cell array if searchmode is 1!');
end

%% Find the indices
if strcmpi(searchmode, 'exact')             % String must be exact

    % Construct a cell array with the same size as cellarray 
    %   but with str duplicated throughout
    str_cell = cell(size(cellarray));    % a cell array to store copies of str
    for k = 1:numel(cellarray)
        str_cell{k} = str;        % make the kth elementhe same as str
    end

    % Find indices that corresponds to str exactly in cellarray, 
    %   case-insensitive if IgnoreCase is set to true
    if ignorecase
        indices = find(cellfun(@strcmpi, cellarray, str_cell));
    else
        indices = find(cellfun(@strcmp, cellarray, str_cell));
    end

elseif strcmpi(searchmode, 'substrings')    % String can be a substring 
                                            % or a cell array of substrings

    % Convert each string to lower case if IgnoreCase is set to true
    if ignorecase
        cellarray = cellfun(@lower, cellarray, 'UniformOutput', false);
    end

    % Find indices that contain str in cellarray
    if iscell(str)        % if str is a cell array of substrings
        % Convert each substring to lower case if IgnoreCase is set to true
        if ignorecase
            str = cellfun(@lower, str, 'UniformOutput', false);
        end

        % Find the indices that contain each substring
        numstrs = numel(str);
        indices_eachstr = cell(1, numstrs);
        for k = 1:numstrs
            indicesarray = strfind(cellarray, str(k));    
            indices_eachstr{k} = find(~cellfun(@isempty, indicesarray));
        end

        % Find the indices that contain all substrings by intersection
        indices = intersectm(indices_eachstr);    
    
    else                    % if str is a single substring
        % Convert substring to lower case if IgnoreCase is set to true
        if ignorecase
            str = lower(str);
        end

        % Find the indices that contain the substring
        indicesarray = strfind(cellarray, str);
        indices = find(~cellfun(@isempty, indicesarray));
    end

end

%% Return the elements too
if ~isempty(indices) && numel(indices) > 1
    elements = cellarray(indices);
elseif ~isempty(indices)
    elements = cellarray{indices};
else
    elements = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

function indices = find_in_strings(str, cellarray, searchmode, ignorecase)

elseif isempty(str) || isempty(cellarray)
    error('First two inputs cannot be empty!');
elseif ~ischar(str) && ~(iscell(str) && min(cellfun(@ischar, str)))
    error('First input must be either a character array or a cell array of character arrays!');
elseif ~iscell(cellarray)
    error('Second input must be a cell array of character arrays!');
elseif nargin >= 3 && ~isempty(searchmode) && ~(searchmode == 1 || searchmode == 2)
    error('Search mode out of range!');
elseif nargin >= 4 && ~isempty(ignorecase) && ~islogical(ignorecase) && ~islogical(ignorecase)

%% Set defaults for optional arguments
if nargin < 3
    searchmode = 2;
end

addParameter(iP, 'SearchMode', 'substrings', ...    % the search mode
    @(x) assert(ischar(x), 'SearchMode must be a either ''exact'' or ''substrings''!'));
addParameter(iP, 'IgnoreCase', false, ...        % whether to ignore any differences in letter case
    @(x) assert(islogical(x) || isequal(x, 1) || isequal(x, 0), 'IgnoreCase must be 1/0 or true/false!'));
addParameter(iP, 'IgnoreCase', false, @islogical);    % whether to ignore any differences in letter case

%}
