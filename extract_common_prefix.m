function prefix = extract_common_prefix (strs, varargin)
%% Extracts the common prefix of a cell array of strings
% Usage: prefix = extract_common_prefix (strs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       extract_common_prefix({'a_b_c', 'a_b_c.m', 'a_b_c_d.m'})
%       extract_common_prefix({'a+b+c', 'a+b+c+d-m'}, 'Delimiter', '+')
%       extract_common_prefix({'a_b_c', 'a_b_c.m'}, 'KeepDelimiter', true)
%       extract_common_prefix({'a_c', 'b_d'})
%       extract_common_prefix({'a_c', 'b_d'}, 'KeepDelimiter', true)
%       extract_common_prefix({'', ''}, 'KeepDelimiter', true)
%       extract_common_prefix('a_b_c')
%       extract_common_prefix('a_b_c', 'KeepDelimiter', true)
%       extract_common_prefix({'a1_c', 'a1_d'}, 'RegExp', '\w')
%
% Outputs:
%       prefix      - the common prefix
%                   specified as a character vector
% Arguments:
%       strs        - strings
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'Delimiter': delimiter used
%                   must be a string scalar or a character vector
%                   default == '_'
%                   - 'KeepDelimiter': whether to keep the trailing delimiter
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SuffixInstead': extract common suffix instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RegExp': regular expression to match
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%                   default == none
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_substrings.m
%       cd/extract_subvectors.m
%       cd/isemptycell.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/extract_common_directory.m
%       cd/extract_common_suffix.m
%       cd/extract_fileparts.m
%       cd/extract_substrings.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/parse_all_abfs.m
%       cd/parse_current_family.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m
%       cd/plot_raster.m
%       cd/plot_raw_multiunit.m
%       cd/plot_spike_density_multiunit.m
%       cd/plot_table.m

% File History:
% 2018-12-26 Moved from extract_common_directory.m
% 2018-12-26 Added 'Delimiter' and 'SuffixInstead' as optional arguments
% 2018-12-27 Added 'KeepDelimiter' as an optional argument
% 2019-08-23 Now uses isemptycell.m
% 2019-12-21 Fixed the case when there is only one string
% 2019-12-21 Added 'RegExp' as an optional argument 
% 

%% Hard-coded parameters

%% Default values for optional arguments
delimiterDefault = '_';
keepDelimiterDefault = false;   % don't keep the preceding delimiter by default
suffixInsteadDefault = false;
regExpDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'strs', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs must be a character array or a string array ', ...
        'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'KeepDelimiter', keepDelimiterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SuffixInstead', suffixInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RegExp', regExpDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['regExp must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));

% Read from the Input Parser
parse(iP, strs, varargin{:});
delimiter = iP.Results.Delimiter;
keepDelimiter = iP.Results.KeepDelimiter;
suffixInstead = iP.Results.SuffixInstead;
regExp = iP.Results.RegExp;

%% Preparation
if suffixInstead
    alignMethod = 'rightadjust';
else
    alignMethod = 'leftadjust';
end
        
%% Do the job
% If empty, return empty
if isemptycell(strs)
    prefix = '';
    return
end

% Only do anything complicated if there is more than one string
if ischar(strs) || isstring(strs) && numel(strs) == 1
    prefix = strs;
elseif numel(strs) == 1
    prefix = strs{1};
else
    % Split all strings by delimiter to get parts
    %   Note: split() can take both a cell and a string as the argument
    %           and returns a column cell array
    parts = arrayfun(@(x) split(x, delimiter), strs, 'UniformOutput', false);

    % Extract the same number of elements from each cell array
    partsAligned = extract_subvectors(parts, 'AlignMethod', alignMethod, ...
                                        'TreatCellStrAsArray', true);

    % Place all parts together in a 2-D cell array
    %   Each column is an original string
    %   Each row is a level
    partsArray = horzcat(partsAligned{:});

    % Separate parts by level
    partsByLevel = extract_rows(partsArray);

    % Find the number of unique elements in each row
    nUniqueEachLevel = count_unique_elements(partsByLevel);

    if suffixInstead
        % Find the last row that has more than one unique element
        levelLastDifference = find(nUniqueEachLevel > 1, 1, 'last');

        % If every row has more than one unique element, 
        %   the common suffix is empty
        if levelLastDifference == numel(nUniqueEachLevel)
            prefix = '';
            return
        end

        % Use the next row number
        if isempty(levelLastDifference) 
            levelFirstCommon = 1;
        else
            levelFirstCommon = levelLastDifference + 1;
        end

        % Construct the common suffix
        tempCell = join(partsAligned{1}(levelFirstCommon:end), delimiter);
        prefix = tempCell{1};
    else
        % Find the first row that has more than one unique element
        levelFirstDifference = find(nUniqueEachLevel > 1, 1, 'first');

        % If every row has more than one unique element, 
        %   the common prefix is empty
        if levelFirstDifference == 1
            prefix = '';
            return
        end

        % Use the previous row number
        if isempty(levelFirstDifference)
            levelLastCommon = numel(nUniqueEachLevel);
        else
            levelLastCommon = levelFirstDifference - 1;
        end

        % Construct the common prefix
        tempCell = join(partsAligned{1}(1:levelLastCommon), delimiter);
        prefix = tempCell{1};
    end
end

% Match regular expression if requested
if ~isempty(regExp)
    prefix = extract_substrings(prefix, 'RegExp', regExp);
end

% Add the delimiter if requested
if keepDelimiter
    if suffixInstead
        prefix = [delimiter, prefix];
    else
        prefix = [prefix, delimiter];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rowsExtracted = extract_rows(cellArray)
%% Extracts rows from a 2D cell array
% TODO: Pull out to its own function
% TODO: Think about how to make a extract_columns work the same way
%       with an optional argument

% Count the number of levels
nRows = size(cellArray, 1);

% Separate parts by level
rowsExtracted = arrayfun(@(x) cellArray(x, :), transpose(1:nRows), ...
                            'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nUnique = count_unique_elements(vecs)
%% Counts the number of unique elements in each vector
% TODO: Pull out to its own function
% TODO: Use extract_subvectors with a 'Unique' option

% Extract unique elements
uniqueElements = cellfun(@unique, vecs, 'UniformOutput', false);

% Count the numbers
nUnique = count_samples(uniqueElements);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

pathParts = cellfun(@(x) split(x, filesep), directories, 'UniformOutput', false);

nParts = cellfun(@numel, pathParts);

% Initialize the last common part to minNParts
ctLastCommon = minNParts;

% Run through all path parts until they become different
for iNPart = 1:minNParts
    % Get this part from all strings
    thisPart = cellfun(@(x) x{iNPart}, pathParts, 'UniformOutput', false);

    % Count the number of unique parts
    nUniqueParts = numel(unique(thisPart));

    % If the number of unique parts is not one, make the previous part the last
    %   common part and exit the loop
    if nUniqueParts ~= 1
        ctLastCommon = iNPart - 1;
        break
    end
end

tempCell = join(pathParts{1}(1:ctLastCommon), filesep);
prefix = tempCell{1};

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
