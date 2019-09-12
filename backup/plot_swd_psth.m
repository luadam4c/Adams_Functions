function handles = plot_swd_psth (varargin)
%% Plots a peri-stimulus time histogram from all gas_pulses.csv files and SWDs.csv files in a directory (unfinished)
% Usage: handles = plot_swd_psth (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       plot_swd_psth('Directory', '/media/shareX/2019octoberR01/Figures/Figure1d')
%
% Outputs:
%       handles     - TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'FirstOnly': whether to take only the first window 
%                                   from each file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': directory to look for SWD and stim table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'RelativeTimeWindowMin': relative time window
%                   must be a 2-element numeric vector
%                   default == interStimInterval * 0.5 * [-1, 1]
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == TODO
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == TODO
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == {'png', 'epsc2'}
%                   - Any other parameter-value pair for plot_psth()
%
% Requires:
%       cd/all_files.m
%       cd/create_label_from_sequence.m
%       cd/extract_distinct_fileparts.m
%       cd/extract_fileparts.m
%       cd/plot_psth.m
%
% Used by:
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2019-09-10 Created by Adam Lu
% TODO: Use load_matching_sheets.m
% 

%% Hard-coded parameters
SEC_PER_MIN = 60;

%   TODO: Create load_gas_pulses.m;
% Note: Must be consistent with parse_iox.m and parse_gas_trace.m
stimTableSuffix = '_gas_pulses';

% Note: Must be consistent with all_swd_sheets.m
swdTableSuffix = '_SWDs';

pathBase = '';
sheetType = 'csv';
figSuffix = '';

%% Default values for optional arguments
firstOnlyDefault = false;       % take all windows by default
directoryDefault = '';          % set later
relativeTimeWindowMinDefault = [];
figTitleDefault = '';           % set later
figNameDefault = '';            % set later
figTypesDefault = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FirstOnly', firstOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RelativeTimeWindowMin', relativeTimeWindowMinDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
firstOnly = iP.Results.FirstOnly;
directory = iP.Results.Directory;
relTimeWindowMin = iP.Results.RelativeTimeWindowMin;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_psth() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default figure suffix
if isempty(figSuffix)
    if firstOnly
        figSuffix = 'first_windows';
    else
        figSuffix = 'all_windows';
    end
end

% Set default figure name
if isempty(figName)
    % Extract the directory base
    dirBase = extract_fileparts(directory, 'dirbase');

    if ~isempty(relTimeWindowMin)
        % Create a sequence for create_label_from_sequence
        relTimeWindowSeq = linspace(relTimeWindowMin(1), ...
                                    relTimeWindowMin(2), 2);

        % Create a figure name
        figName = fullfile(directory, [dirBase, '_', ...
                            create_label_from_sequence(relTimeWindowSeq), ...
                            '_', figSuffix]);
    else
        % Create a figure name
        figName = fullfile(directory, [dirBase, '_', figSuffix]);
    end
end

% Set default figure title
if isempty(figTitle)
    if firstOnly
        figTitle = '# of SWDs around first stimulus starts';
    else
        figTitle = '# of SWDs around all stimulus starts';
    end
end

%% Do the job
% Get all stim pulse files
[~, stimPaths] = ...
    all_files('Directory', directory, 'Keyword', pathBase, ...
                'Suffix', stimTableSuffix, 'Extension', sheetType);

% Extract distinct prefixes
distinctPrefixes = extract_distinct_fileparts(stimPaths);

% Look for corresponding SWD spreadsheet files
[~, swdPaths] = ...
    cellfun(@(x) all_files('Prefix', x, 'Directory', directory, ...
                    'Suffix', swdTableSuffix, 'Extension', sheetType), ...
            distinctPrefixes, 'UniformOutput', false);

% Read all tables
[stimTables, swdTables] = ...
    argfun(@(x) cellfun(@readtable, x, 'UniformOutput', false), ...
            stimPaths, swdPaths);

% Extract all start times in seconds
[stimStartTimesSec, swdStartTimesSec] = ...
    argfun(@(x) cellfun(@(y) y.startTime, x, 'UniformOutput', false), ...
            stimTables, swdTables);

% Restrict to the first events only if requested
if firstOnly
    % Note: TODO: Use extract_elements.m with 'UniformOutput', false
    stimStartTimesSec = cellfun(@(x) x(1), stimStartTimesSec, ...
                                'UniformOutput', false);
end

% Convert to minutes
[stimStartTimesMin, swdStartTimesMin] = ...
    argfun(@(x) cellfun(@(y) y / SEC_PER_MIN, x, 'UniformOutput', false), ...
            stimStartTimesSec, swdStartTimesSec);

%% Plot the peri-stimulus time histogram
handles = plot_psth('EventTimes', swdStartTimesMin, ...
                    'StimTimes', stimStartTimesMin, ...
                    'XLabel', 'Time (min)', ...
                    'RelativeTimeWindow', relTimeWindowMin, ...
                    'FigTitle', figTitle, ...
                    'FigName', figName, 'FigTypes', figTypes, ...
                    otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
