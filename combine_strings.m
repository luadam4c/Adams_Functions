function finalString = combine_strings (varargin)
%% Constructs a final string based on optional substrings and/or Name-Value pairs
% Usage: finalString = combine_strings (strs (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       combine_strings
%       combine_strings('Substrings', {'_yes__', '_no_'})
%       combine_strings('Substrings', {'_yes__', '_no_'}, 'ForceClean', false)
%       combine_strings('Substrings', {'', 'test', ''})
%       combine_strings('Substrings', {{'funny', 'boy'}, 'test'})
%       combine_strings('Substrings', {{'funny', 'boy'}, {'high', 'low'}})
%       combine_strings('NameValuePairs', {{'a', 'b'}, [1, 2]})
%       combine_strings('Substrings', {'yes', 'no'}, 'NameValuePairs', {{'a', 'b'}, [1, 2]})
%
% Outputs:
%       finalString    - a string (may be empty) that is a final substring
%
% Arguments:
%       strs     - (opt) either 'Substrings' or 'NameValuePairs'
%       varargin    - 'Delimiter': delimiter used
%                   must be a string scalar or a character vector
%                   default == '_'
%                   - 'ForceClean': whether the delimiter is 
%                                   not to be repeated between substrings
%                                   and trimmed at either ends
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'BeginWithDelimiter': whether to begin with the delimiter
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'EndWithDelimiter': whether to end with the delimiter
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Substrings': substring(s) to combine
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                       or a cell array of those
%                   default == ''
%                   - 'NameValuePairs': Name-Value pairs that are changed
%                   must be a 2-element cell array whose first element 
%                       is a string/char array or cell array 
%                       and whose second element is a numeric array
%                   default == {'', NaN}
%        
% Requires:
%       cd/isemptycell.m
%       cd/force_string_end.m
%       cd/force_string_start.m
%
% Used by:
%       cd/construct_fullpath.m
%       cd/m3ha_autocorrelogram.m
%       cd/m3ha_network_raster_plot.m
%       cd/m3ha_neuron_choose_best_params.m
%       /media/adamX/RTCl/raster_plot.m

% File History:
% 2017-05-04 Moved from construct_fullfilename.m
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2019-11-28 Now accepts a cell array of strings for parts
% 2019-11-28 Added 'ForceClean' as an optional argument 
%           (default == true), where the delimiter is checked not to be repeated
% TODO Modify the optional argument 'strs' that could be 
%               interpreted as either 'Substrings' or 'NameValuePairs'
% TODO: Change specification of NameValuePairs to just one cell array
%       or a structure and use struct2arglist.m

%% Hard-coded parameters

%% Default values for optional arguments
strsDefault = {};
delimiterDefault = '_';
forceCleanDefault = true;
beginWithDelimiterDefault = false;
endWithDelimiterDefault = false;
substringsDefault = '';
nameValuePairsDefault = {'', NaN};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add optional arguments to the input Parser
addOptional(iP, 'strs', strsDefault, ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ForceClean', forceCleanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'BeginWithDelimiter', beginWithDelimiterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'EndWithDelimiter', endWithDelimiterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Substrings', substringsDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x) || iscell(x), ...
        ['Substrings must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'NameValuePairs', nameValuePairsDefault, ...
    @(x) assert(iscell(x) && numel(x) == 2, ...
                'NameValuePairs must be a 2-element cell array!'));

% Read from the input Parser
parse(iP, varargin{:});
strs = iP.Results.strs;
delimiter = iP.Results.Delimiter;
forceClean = iP.Results.ForceClean;
beginWithDelimiter = iP.Results.BeginWithDelimiter;
endWithDelimiter = iP.Results.EndWithDelimiter;
substrings = iP.Results.Substrings;
nameValuePairs = iP.Results.NameValuePairs;

%% Determine what the first argument is
if ~isempty(strs)
    substrings = strs;
    % if iscellstr(strs) || TODO
    %     substrings = strs;
    % else
    %     nameValuePairs = strs;
    % end
end

%% Find all substrings
if isempty(substrings) && isempty(nameValuePairs{1})
    allSubstrings = '';
else
    % Initialize a cell array for all substrings
    numSuffixes = 0;
    if iscell(substrings)
        numSuffixes = numSuffixes + numel(substrings);
    elseif ~isempty(substrings)
        numSuffixes = numSuffixes + 1;
    end
    if iscell(nameValuePairs{1})
        numSuffixes = numSuffixes + numel(nameValuePairs{1});
    elseif ~isempty(nameValuePairs{1})
        numSuffixes = numSuffixes + 1;
    end
    allSubstrings = cell(1, numSuffixes);             % stores all substrings
    ct = 0;

    % Add premade substrings if premade substrings are provided
    if ~isempty(substrings)
        if iscell(substrings)
            for s = 1:numel(substrings)
                ct = ct + 1;
                allSubstrings{ct} = substrings{s};
            end
        else
            ct = ct + 1;
            allSubstrings{ct} = substrings;
        end
    end

    % Add premade Name-Value pairs if 
    %   Name-Value pairs that are changed are provided
    if ~isempty(nameValuePairs{1})
        if iscell(nameValuePairs{1})
            % If there might be more than one Name-Value pairs provided,
            %   add iteratively
            for p = 1:numel(nameValuePairs{1})
                ct = ct + 1;
                allSubstrings{ct} = [nameValuePairs{1}{p}, delimiter, ...
                                    num2str(nameValuePairs{2}(p))];
            end
        else
            % If there is only one Name-Value pair provided,
            %   add the pair only
            ct = ct + 1;
            allSubstrings{ct} = [nameValuePairs{1}, delimiter, ...
                                num2str(nameValuePairs{2})];
        end
    end
end

%% Construct final substring
if all(isemptycell(allSubstrings))            % if nothing provided
    % Final substring is empty too
    finalString = '';
else
    if iscell(allSubstrings)
        % If there might be more than one substrings provided,
        %   construct final substring by concatenating all substrings 
        %   together with delimiter
        finalString = strjoin_custom(allSubstrings, delimiter, forceClean);
    else
        % If there is only one substring provided, 
        %   the final substring is this substring
        finalString = allSubstrings;
    end
end

%% Force the substring to start with delimiter if requested
if beginWithDelimiter
    finalString = force_string_start(finalString, delimiter, ...
                                    'OnlyIfNonempty', true);
end

%% Force the substring to end with delimiter if requested
if endWithDelimiter
    finalString = force_string_end(finalString, delimiter, ...
                                    'OnlyIfNonempty', true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function joinedStrs = strjoin_custom(subStrs, delimiter, forceClean)
%% Joins strings 

% Modify parts
if forceClean
    % Strip the leading and trailing delimiters
    subStrs = strip_custom(subStrs, delimiter);

    % Force each nonempty substr to end with one delimiter
    subStrsWithDelimiter = force_string_end(subStrs, delimiter, ...
                                            'OnlyIfNonempty', true);
else
    % Attach delimiter to the end of each substr
    subStrsWithDelimiter = strcat(subStrs, delimiter);
end

% Concatenate all substrings
if iscell(subStrsWithDelimiter)
    joinedStrsWithEndDelimiter = strcat(subStrsWithDelimiter{:});
else
    joinedStrsWithEndDelimiter = strcat(subStrsWithDelimiter(:));
end

% Remove the ending delimiter
joinedStrs = remove_if_end(joinedStrsWithEndDelimiter, delimiter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outStr = strip_custom (inStr, stripCharacter)
%% Strips a character from all strings in first argument
% Example: 
%           strip_custom({{'a1', 'b1'}, 'c1'}, '1')
% TODO: Pull out as a strip_custom.m

if iscell(inStr)
    outStr = cellfun(@(x) strip_custom(x, stripCharacter), inStr, ...
                        'UniformOutput', false);
else
    outStr = strip(inStr, stripCharacter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outStr = remove_if_end (inStr, subStr)
%% Removes a substring from strings if it is at the end
% TODO: Pull out as a function remove_if_end.m
% TODO: make remove_if_start.m

if iscell(inStr)
    % Do this for each string separately
    outStr = cellfun(@(x) remove_if_end(x, subStr), ...
                        inStr, 'UniformOutput', false);
elseif isstring(inStr) && ...
        numel(inStr) > 1
    % Do this for each string separately
    outStr = arrayfun(@(x) remove_if_end(x, subStr), ...
                        inStr, 'UniformOutput', true);
else
    % Compute the number of characters in each string
    nCharsInStr = strlength(inStr);
    nCharsSubStr = strlength(subStr);

    % Only take out subStr if it exists
    if nCharsInStr >= nCharsSubStr
        % Count the number of characters to extract
        nCharsToExtract = nCharsInStr - nCharsSubStr;

        % Extract everything before the next character
        outStr = extractBefore(inStr, nCharsToExtract + 1);
    else
        % Don't extract anything
        outStr = inStr;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

addParameter(iP, 'NameValuePairs', nameValuePairsDefault, ...
    @(x) assert(iscell(x) && numel(x) == 2 ...
            && (ischar(x{1}) || iscell(x{1}) || isstring(x{1})) ...
            && isnumeric(x{2}), ...
        ['NameValuePairs must be a 2-element cell array whose ', ...
            'first element is a string/char array or cell array ', ...
            'and whose second element is a numeric array!']));

% If there might be more than one substrings provided,
%   construct final substring by concatenating all substrings 
%   together with '_'
finalString = allSubstrings{1};
if numel(allSubstrings) > 1
    for s = 2:numel(allSubstrings)
        finalString = [finalString, '_', allSubstrings{s}];
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
