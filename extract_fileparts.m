function parts = extract_fileparts (paths, partType, varargin)
%% Extracts directories, bases, extensions, distinct parts or the common directory from file paths, treating any path without an extension as a directory
% Usage: parts = extract_fileparts (paths, partType, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [~, paths] = all_files('Directory', pwd);
%       paths = {'a/b/c.m', 'a/b/c'}
%       extract_fileparts(paths, 'directory')
%       extract_fileparts(paths, 'dirbase')
%       extract_fileparts(paths, 'parentdir')
%       extract_fileparts(paths, 'filename')
%       extract_fileparts(paths, 'filebase')
%       extract_fileparts(paths, 'name')
%       extract_fileparts(paths, 'base')
%       extract_fileparts(paths, 'pathbase')
%       extract_fileparts(paths, 'ext')
%       extract_fileparts(paths, 'commondirectory')
%       extract_fileparts(paths, 'commonprefix')
%       extract_fileparts(paths, 'commonsuffix')
%       extract_fileparts(paths, 'distinct')
%       sliceNames = extract_fileparts(paths, 'base', 'RegExp', '.*slice[0-9]*');
%
% Outputs:
%       parts       - parts extracted
%                   specified as a character array 
%                       or a cell array of character arrays
%
% Arguments:
%       paths       - file paths (absolute or relative)
%                   Note: this must contain the extension if it is a file
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       partType    - type of the file part to extract
%                   must be an unambiguous, case-insensitive match to one of:
%                       'commondirectory' - common directory across file(s)
%                       'commonprefix'    - common prefix across file(s)
%                       'commonsuffix'    - common suffix across file(s)
%                       'distinct'  - distinct parts across file(s)
%                       'directory' - directory containing the file(s)
%                       'dirbase'   - directory base containing the file(s)
%                       'parentdir' - directory one level up
%                       'name'      - file name with the extension
%                       'base'      - file base name without the extension
%                       'pathbase'  - full file path without the extension
%                       'extension' - file extension including the leading '.'
%       varargin    - 'Delimiter': delimiter used for file suffixes
%                   must be a string scalar or a character vector
%                   default == '_'
%                   - 'RegExp': regular expression to match
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%                   default == none
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_new_mscript.m
%       cd/extract_common_directory.m
%       cd/extract_common_prefix.m
%       cd/extract_common_suffix.m
%       cd/extract_distinct_fileparts.m
%
% Used by:
%       cd/all_dependent_functions.m
%       cd/all_file_bases.m
%       cd/archive_dependent_scripts.m
%       cd/atfwrite.m
%       cd/combine_abf_data.m
%       cd/combine_data_from_same_slice.m
%       cd/combine_swd_sheets.m
%       cd/combine_variables_across_tables.m
%       cd/copy_into.m
%       cd/create_pleth_EEG_movies.m
%       cd/decide_on_filebases.m
%       cd/extract_common_directory.m
%       cd/find_matching_files.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_network_single_neuron.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_plot_violin.m
%       cd/m3ha_select_sweeps.m
%       cd/parse_atf_swd.m
%       cd/parse_current_family.m
%       cd/parse_iox.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_measures.m
%       cd/plot_relative_events.m
%       cd/plot_small_chevrons.m
%       cd/plot_swd_histogram.m
%       cd/plot_table.m
%       cd/plot_traces.m
%       cd/plot_traces_spike2_mat.m
%       cd/save_all_zooms.m
%       cd/update_file_base_in_matfiles.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2018-12-18 Created by Adam Lu
% 2018-12-26 Added 'commonsuffix' as a part type
% 2018-12-27 Moved code to extract_distinct_fileparts.m
% 2019-03-14 Added 'commonprefix' as a part type
% 2019-04-07 Added 'dirbase' as a part type
% 2019-04-07 Added 'RegExp' as an optional argument
% 2019-04-07 Added 'name' as a part type
% 2019-08-12 Added 'pathbase' as a part type
% 2019-09-10 Fixed 'dirbase' when the input is a cell array
% 2019-09-30 Added 'parentdir' as a part type
% 2019-12-22 Now makes 'base' different from 'filebase'
%                       and 'name' different from 'filename'
% TODO: Make the first argument accept a files structure array too
% 

%% Hard-coded parameters
validPartTypes = {'commondirectory', 'commonprefix', 'commonsuffix', ...
                    'distinct', 'directory', 'dirbase', 'parentdir', ...
                    'name', 'base', 'filename', 'filebase', ...
                    'pathbase', 'extension'};

%% Default values for optional arguments
delimiterDefault = '_';
regExpDefault = '';

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
addRequired(iP, 'paths', ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['paths must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'partType', ...
    @(x) any(validatestring(x, validPartTypes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RegExp', regExpDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['regExp must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));

% Read from the Input Parser
parse(iP, paths, partType, varargin{:});
delimiter = iP.Results.Delimiter;
regExp = iP.Results.RegExp;

% Validate partType
partType = validatestring(partType, validPartTypes);

%% Preparation


%% Do the job
switch partType
case {'directory', 'dirbase', 'parentdir', ...
        'name', 'base', 'filename', 'filebase', 'pathbase', 'extension'}
    parts = extract_simple_fileparts(paths, partType);
case 'commondirectory'
    parts = extract_common_directory(paths, varargin{:});
case {'commonprefix', 'commonsuffix'}
    % First, extract file bases
    fileBases = extract_simple_fileparts(paths, 'base');

    % Next, extract file prefixes or suffixes
    switch partType
        case 'commonprefix'
            parts = extract_common_prefix(fileBases, 'Delimiter', delimiter);
        case 'commonsuffix'
            parts = extract_common_suffix(fileBases, 'Delimiter', delimiter);
    end
case 'distinct'
    parts = extract_distinct_fileparts(paths, 'Delimiter', delimiter);
otherwise
    error('partType unrecognized!!');
end

% Match regular expression if provided
if ~isempty(regExp)
    parts = regexp(parts, regExp, 'match', ...
                    'once', 'ignorecase', 'emptymatch', 'dotexceptnewline');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fileDir, fileBase, fileExtension] = ...
                fileparts_custom (filePath, useOriginal)
%% Same as fileparts but treats anything without an extension as a directory

if useOriginal
    [fileDir, fileBase, fileExtension] = fileparts(filePath);
else
    % Use the original fileparts to get the extension
    [fileDirTentative, fileBaseTentative, fileExtension] = fileparts(filePath);

    % Check if there is an extension
    if isempty(fileExtension)
        % If there is no extension, it is a directory
        fileDir = fullfile(fileDirTentative, fileBaseTentative);
        fileBase = '';
    else
        % If there is an extension, it is a file
        fileDir = fileDirTentative;
        fileBase = fileBaseTentative;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parts = extract_simple_fileparts(paths, partType)
% Separate the paths with fileparts_custom()

% Decide whether to use the original fileparts directly
switch partType
    case {'base', 'name', 'pathbase'}
        useOriginal = true;
    otherwise
        useOriginal = false;
end

% Separate file parts
if iscell(paths)
    [fileDirs, fileBases, fileExtensions] = ...
        cellfun(@(x) fileparts_custom(x, useOriginal), ...
                paths, 'UniformOutput', false);
else
    [fileDirs, fileBases, fileExtensions] = ...
        fileparts_custom(paths, useOriginal);
end

% Return results
switch partType
    case 'directory'
        parts = fileDirs;
    case 'dirbase'
        parts = extract_simple_fileparts(strcat(fileDirs, '.dum'), 'base');
    case 'parentdir'
        parts = extract_simple_fileparts(strcat(fileDirs, '.dum'), 'directory');
    case {'filename', 'name'}
        parts = strcat(fileBases, fileExtensions);
    case {'filebase', 'base'}
        parts = fileBases;
    case 'pathbase'
        parts = fullfile(fileDirs, fileBases);        
    case 'extension'
        parts = fileExtensions;
    otherwise
        error('partType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
