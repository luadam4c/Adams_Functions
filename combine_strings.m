function finalString = combine_strings (varargin)
%% Constructs a final string based on optional substrings and/or Name-Value pairs
% Usage: finalString = combine_strings (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       combine_strings
%       combine_strings('Substrings', {'_yes__', '_no_'})
%       combine_strings('Substrings', {'_yes__', '_no_'}, 'ForceClean', false)
%       combine_strings('Substrings', {{'funny', 'boy'}, 'test'})
%       combine_strings('Substrings', {{'funny', 'boy'}, {'high', 'low'}})
%       combine_strings('NameValuePairs', {{'a', 'b'}, [1, 2]})
%       combine_strings('Substrings', {'yes', 'no'}, 'NameValuePairs', {{'a', 'b'}, [1, 2]})
%
% Outputs:
%       finalString    - a string (may be empty) that is a final substring
%
% Arguments:
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
% TODO: Change specification of NameValuePairs to just one cell array
%       or a structure and use struct2arglist.m
% 2019-11-28 Added 'ForceClean' as an optional argument 
%           (default == true), where the delimiter is checked not to be repeated

%% Hard-coded parameters

%% Default values for optional arguments
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
delimiter = iP.Results.Delimiter;
forceClean = iP.Results.ForceClean;
beginWithDelimiter = iP.Results.BeginWithDelimiter;
endWithDelimiter = iP.Results.EndWithDelimiter;
substrings = iP.Results.Substrings;
nameValuePairs = iP.Results.NameValuePairs;

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
if isempty(allSubstrings)            % if nothing provided
    % Final substring is empty too
    finalString = allSubstrings;
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
nCharsToExtract = strlength(joinedStrsWithEndDelimiter) - strlength(delimiter);
joinedStrs = extractBefore(joinedStrsWithEndDelimiter, nCharsToExtract + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subStrs = strip_custom(subStrs, delimiter)

if iscell(subStrs)
    subStrs = cellfun(@(x) strip_custom(x, delimiter), subStrs, ...
                        'UniformOutput', false);
else
    subStrs = strip(subStrs, delimiter);
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
