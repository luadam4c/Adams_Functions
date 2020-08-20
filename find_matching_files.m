function varargout = find_matching_files (fileStrs, varargin)
%% Finds matching files from file strings
% Usage: [files, fullPaths, distinctParts] = find_matching_files (fileStrs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [~, matPaths] = all_files('Ext', 'mat');
%       [csvFiles, csvPaths] = find_matching_files(matPaths, 'Extension', 'csv');
%       [wmvFiles, wmvPaths] = find_matching_files(matPaths, 'Extension', 'wmv');
%
% Outputs:
%       files       - file structure(s) for the files
%                   specified as a structure array with fields:
%                       name
%                       folder
%                       date
%                       bytes
%                       isdir
%                       datenum
%       fullPaths   - full path(s) to the files
%                   specified as a column cell array of character vectors
%       distinctParts   - distinct parts between different files
%                   specified as a column cell array of character vectors
%
% Arguments:
%       fileStrs   - file strings to match (can be full paths)
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'PartType': part type to match
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'prefix'    - match the prefix
%                       'keyword'   - match any part of the file name
%                       'suffix'    - match the suffix
%                   default == 'keyword'
%                   - 'ExtractDistinct': whether to extract distinct parts
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ReturnEmpty': whether to return an empty string 
%                                   if not matched
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for all_files()
%
% Requires:
%       cd/all_files.m
%       cd/create_error_for_nargin.m
%       cd/extract_distinct_fileparts.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/find_first_match.m
%       cd/force_string_end.m
%
% Used by:
%       cd/combine_swd_resp_data.m
%       cd/create_pleth_EEG_movies.m
%       cd/create_power_tables.m
%       cd/read_matching_sheets.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_simulate_population.m
%       cd/plot_relative_events.m
%       cd/plot_traces_spike2_mat.m

% File History:
% 2019-09-25 Created by Adam Lu
% 2019-09-30 Now maintains character vectors as character vectors
% 2019-10-15 Added 'ForceCellOutput' as an optional argument
% 2019-12-20 Changed default extractDistinct to false
% 2020-02-02 Added 'ReturnEmpty' as an optional argument
% 2020-02-09 Fixed bug if file not found and returning empty
% TODO: Add 'Delimiter' as an optional argument
% TODO: 'MaxNum' not always 1
% 

%% Hard-coded parameters
validPartTypes = {'prefix', 'keyword', 'suffix'};

%% Default values for optional arguments
partTypeDefault = 'keyword';
extractDistinctDefault = false; % don't extract distinct parts by default
returnEmptyDefault = false;
forceCellOutputDefault = false; % don't force output as a cell array by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'fileStrs', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PartType', partTypeDefault, ...
    @(x) any(validatestring(x, validPartTypes)));
addParameter(iP, 'ExtractDistinct', extractDistinctDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ReturnEmpty', returnEmptyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, fileStrs, varargin{:});
partType = validatestring(iP.Results.PartType, validPartTypes);
extractDistinct = iP.Results.ExtractDistinct;
returnEmpty = iP.Results.ReturnEmpty;
forceCellOutput = iP.Results.ForceCellOutput;

% Keep unmatched arguments for the all_files() function
otherArguments = iP.Unmatched;

%% Extract distinct parts
% Force as a cell array
if ischar(fileStrs)
    fileStrs = force_column_cell(fileStrs);
    wasChar = true;
else
    wasChar = false;
end

% Extract distinct file strings
if extractDistinct
    distinctParts = extract_distinct_fileparts(fileStrs);
else
    distinctParts = fileStrs;
end

% Extract the base name within distinct parts
distinctPartsBase = extract_fileparts(distinctParts, 'base');

% Extract the directory name within distinct parts
distinctPartsDir = extractBefore(distinctParts, distinctPartsBase);

% Decide if all files are from the same directory
if ischar(distinctPartsDir)
    commonDir = distinctPartsDir;
    findFilesOneByOne = false;
else
    uniqueDirs = unique(distinctPartsDir);
    if numel(uniqueDirs) == 1
        commonDir = uniqueDirs{1};
        findFilesOneByOne = false;
    else
        findFilesOneByOne = true;
    end
end

%% Do the job
if findFilesOneByOne
    % Find one matching file for each file string
    [filesCell, fullPaths] = ...
        cellfun(@(x, y) all_files('Directory', x, partType, y, 'MaxNum', 1, ...
                            'ForceCellOutput', false, otherArguments), ...
                distinctPartsDir, distinctPartsBase, 'UniformOutput', false);

    % Try to convert to an array
    %   Note: this fails if a cell is empty
    try
        files = cellfun(@(x) x, filesCell);
    catch
        disp([mfilename, ': Some files were not found!']);
        files = filesCell;
    end
else
    % Find all files under the common directory
    [files, fullPaths] = ...
        all_files('Directory', commonDir, 'ForceCellOutput', true, ...
                    otherArguments);

    % Extract just the file bases
    fileBases = extract_fileparts(fullPaths, 'base');

    % Find the first matching file
    indMatched = find_first_match(distinctPartsBase, fileBases, ...
                                        'MatchMode', partType);

    % Extract the corresponding files
    if ~any(isnan(indMatched))
        files = files(indMatched);
        fullPaths = fullPaths(indMatched);
    elseif ~returnEmpty
        error('Some files not matched!');
    else
        files = extract_subvectors(files, 'Indices', indMatched);
        fullPaths = extract_subvectors(fullPaths, 'Indices', indMatched);
    end
end

% Extract the character array if it was one
if wasChar && ~forceCellOutput
    if isempty(fullPaths)
        fullPaths = '';
    else
        fullPaths = fullPaths{1};
    end
end

% Get first output
varargout{1} = files;
varargout{2} = fullPaths;
varargout{3} = distinctParts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
