function finalSuffix = construct_suffix (varargin)
%% Constructs final suffix based on optional suffixes and/or Name-Value pairs
% Usage: finalSuffix = construct_suffix (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       construct_suffix
%       construct_suffix('Suffixes', {'yes', 'no'})
%       construct_suffix('NameValuePairs', {{'a', 'b'}, [1, 2]})
%       construct_suffix('Suffixes', {'yes', 'no'}, 'NameValuePairs', {{'a', 'b'}, [1, 2]})
%
% Outputs:
%       finalSuffix    - a string (may be empty) that is a final suffix
%
% Arguments:
%       varargin    - 'BeginWithDelimiter': whether to begin with '_'
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Suffixes': suffix(ces) to add to filebase
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%                   default == none
%                   - 'NameValuePairs': Name-Value pairs that are changed
%                   must be a 2-element cell array whose first element 
%                       is a string/char array or cell array 
%                       and whose second element is a numeric array
%                   default == none
%        
%
% Used by:
%       cd/construct_fullpath.m
%       cd/force_string_start.m
%       /media/adamX/RTCl/raster_plot.m

% File History:
% 2017-05-04 Moved from construct_fullfilename.m
% 2018-05-08 Changed tabs to spaces and limited width to 80
% TODO: Change specification of NameValuePairs to just one cell array
%       or a structure and use struct2arglist.m

%% Default values for optional arguments
beginWithDelimiterDefault = false;
suffixesDefault = '';
nameValuePairsDefault = {'', NaN};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the input Parser
addParameter(iP, 'BeginWithDelimiter', beginWithDelimiterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Suffixes', suffixesDefault, ...
    @(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x), ...
                ['Suffixes must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));
addParameter(iP, 'NameValuePairs', nameValuePairsDefault, ...
    @(x) assert(iscell(x) && numel(x) == 2, ...
                'NameValuePairs must be a 2-element cell array!'));

% Read from the input Parser
parse(iP, varargin{:});
beginWithDelimiter = iP.Results.BeginWithDelimiter;
suffixes = iP.Results.Suffixes;
nameValuePairs = iP.Results.NameValuePairs;

%% Find all suffixes
if isempty(suffixes) && isempty(nameValuePairs{1})
    allsuffixes = '';
else
    % Initialize a cell array for all suffixes
    numSuffixes = 0;
    if iscell(suffixes)
        numSuffixes = numSuffixes + numel(suffixes);
    elseif ~isempty(suffixes)
        numSuffixes = numSuffixes + 1;
    end
    if iscell(nameValuePairs{1})
        numSuffixes = numSuffixes + numel(nameValuePairs{1});
    elseif ~isempty(nameValuePairs{1})
        numSuffixes = numSuffixes + 1;
    end
    allsuffixes = cell(1, numSuffixes);             % stores all suffixes
    ct = 0;

    % Add premade suffixes if premade suffixes are provided
    if ~isempty(suffixes)
        if iscell(suffixes)
            for s = 1:numel(suffixes)
                ct = ct + 1;
                allsuffixes{ct} = suffixes{s};
            end
        else
            ct = ct + 1;
            allsuffixes{ct} = suffixes;
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
                allsuffixes{ct} = [nameValuePairs{1}{p}, '_', ...
                                    num2str(nameValuePairs{2}(p))];
            end
        else
            % If there is only one Name-Value pair provided,
            %   add the pair only
            ct = ct + 1;
            allsuffixes{ct} = [nameValuePairs{1}, '_', ...
                                num2str(nameValuePairs{2})];
        end
    end
end

%% Construct final suffix
if isempty(allsuffixes)            % if nothing provided
    % Final suffix is empty too
    finalSuffix = allsuffixes;
elseif ~isempty(allsuffixes)        % suffix(ces) is(are) provided
    if iscell(allsuffixes)
        % If there might be more than one suffixes provided,
        %   construct final suffix by concatenating all suffixes 
        %   together with '_'
        finalSuffix = strjoin(allsuffixes, '_');
    else
        % If there is only one suffix provided, 
        %   the final suffix is this suffix
        finalSuffix = allsuffixes;
    end
end

%% Force the suffix to start with '_' if requested
if beginWithDelimiter
    finalSuffix = force_string_start(finalSuffix, '_', 'OnlyIfNonempty', true);
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

% If there might be more than one suffixes provided,
%   construct final suffix by concatenating all suffixes 
%   together with '_'
finalSuffix = allsuffixes{1};
if numel(allsuffixes) > 1
    for s = 2:numel(allsuffixes)
        finalSuffix = [finalSuffix, '_', allsuffixes{s}];
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
