function finalString = combine_strings (varargin)
%% Constructs a final string based on optional substrings and/or Name-Value pairs
% Usage: finalString = combine_strings (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       combine_strings
%       combine_strings('Substrings', {'yes', 'no'})
%       combine_strings('NameValuePairs', {{'a', 'b'}, [1, 2]})
%       combine_strings('Substrings', {'yes', 'no'}, 'NameValuePairs', {{'a', 'b'}, [1, 2]})
%
% Outputs:
%       finalString    - a string (may be empty) that is a final substring
%
% Arguments:
%       varargin    - 'BeginWithDelimiter': whether to begin with the delimiter
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'EndWithDelimiter': whether to end with the delimiter
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Substrings': substring(s) to combine
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
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
%       /media/adamX/RTCl/raster_plot.m

% File History:
% 2017-05-04 Moved from construct_fullfilename.m
% 2018-05-08 Changed tabs to spaces and limited width to 80
% TODO: Change specification of NameValuePairs to just one cell array
%       or a structure and use struct2arglist.m
% TODO: Add 'CheckDelimiter' as an optional argument (default == true), 
%           where the delimiter is checked not to be repeated

%% Hard-coded parameters
% TODO: Make optional arguments
delimiter = '_';

%% Default values for optional arguments
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
addParameter(iP, 'BeginWithDelimiter', beginWithDelimiterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'EndWithDelimiter', endWithDelimiterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Substrings', substringsDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['Substrings must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'NameValuePairs', nameValuePairsDefault, ...
    @(x) assert(iscell(x) && numel(x) == 2, ...
                'NameValuePairs must be a 2-element cell array!'));

% Read from the input Parser
parse(iP, varargin{:});
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
        finalString = strjoin(allSubstrings, delimiter);
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
