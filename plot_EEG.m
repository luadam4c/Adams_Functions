function plot_EEG(abfFileName, varargin)
%% Plots EEG traces from a .abf file and 
% Usage: plot_EEG(abfFileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       abfFileName - file name of the abf file
%                       could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SwdManualFile': TODO: Description of swdManualFileName
%                   must be a TODO
%                   default == TODO
%                   
% Requires:
%       cd/all_files.m
%       cd/parse_abf.m
%       cd/parse_swd_manual.m
%       cd/plot_traces_abf.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
%   2016-09-XX - Created
%   2018-11-21 - Renamed spikewavedetection -> plot_EEG
%   2018-11-21 - Updated to use plot_traces_abf.m
% 

%% Default values for optional arguments
verboseDefault = true;
swdManualFileDefault = '';  % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'abfFileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SwdManualFile', swdManualFileDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, abfFileName, varargin{:});
verbose = iP.Results.Verbose;
swdManualFile = iP.Results.SwdManualFile;

%% Preparation
% Parse the abf file
[parsedParams, parsedData] = ...
    parse_abf(abfFileName, 'Verbose', verbose, 'ExpMode', 'EEG');

% Extract parsed results if there is anything parsed
if isempty(parsedParams)
    return
else
    % Get the file path
    abfPath = parsedParams.abfFullFileName;
end

% Extract the file base
[fileDir, fileBase, ~] = fileparts(abfPath);

% Look for manual SWD files if not provided
if isempty(swdManualFile)
    % Try to look for a .atf file containing the file base in the same directory
    [~, swdManualFullPaths] = all_files('Directory', fileDir, ...
                                        'Keyword', fileBase, ...
                                        'Extension', '.atf');

    % If not successful, try to look for a _Manual_SWDs.csv file containing 
    %   the file base in the same directory
    if isempty(swdManualFullPaths)
        [~, swdManualFullPaths] = all_files('Directory', fileDir, ...
                                            'Keyword', fileBase, ...
                                            'Suffix', '_Manual_SWDs', ...
                                            'Extension', '.csv');
    end

    % If there is more than one file, one choose the first one
    if isempty(swdManualFullPaths)
        swdManualFile = '';
    elseif iscell(swdManualFullPaths)
        swdManualFile = swdManualFullPaths{1};
    else
        swdManualFile = swdManualFullPaths;
    end
end

%% Do the job
% Plot full traces if not already done so
plot_traces_abf(abfPath, 'Verbose', verbose, 'OverWrite', false, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData);

% Plot SWD regions if an SWD file(s) exists
if isfile(swdManualFile)
    % Parse the SWDs from the manual results
    SWDs = parse_swd_manual(swdManualFile);

    % Plot each detected SWD separately
    plot_traces_abf(SWDs.abfPath, 'Verbose', verbose, 'OverWrite', false, ...
            'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
            'TimeStart', SWDs.startTimes, 'TimeEnd', SWDs.endTimes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[alldata, si] = plot_traces_abf_EEG(strcat(cd, '/1_10_2014_cage2_chunk.abf'), 35, 38);
plot_traces_abf_EEG(strcat(cd, '/1_10_2014_cage2_chunk.abf'), 779, 781);
plot_traces_abf_EEG(strcat(cd, '/1_10_2014_cage2_chunk.abf'), 77, 80);
plot_traces_abf_EEG(strcat(cd, '/1_10_2014_cage2_chunk.abf'), 186, 189);

% Load raw data
m = matfile('DBA01_31_2016_cage1.mat');
timeVec = m.t;
channel1 = m.vec1;
channel2 = m.vec2;
channel3 = m.vec3;

table = xlsread('DBA01_31_2016_cage1_Katie.xlsx');
nRows = size(table, 1);     % # of detected SWDs

fprintf(['Check: Duration recorded is %g; 
        duration calculated from start and end times is %g\n'], ...
        duration, endTime - startTime);

% Get the current start time
startTime = table(k, 2) * 60;    % in seconds
% Get the current end time
endTime = table(k, 3) * 60;    % in seconds

duration = table(k, 4);        % in seconds

parfor k = 1:nRows
end

% Count the number of detected SWDs
nRows = height(swdManualTableOfInterest);

% Construct file names for results 
if isempty(swdManualFile)
    swdManualFile = replace(abfPath, '.abf', '_Manual_SWDs.csv');
end

% Construct full path to abf file
[abfPath, pathExists] = construct_and_check_abfpath(abfFileName);
if ~pathExists
    return
end

swdManualFile = replace(abfPath, '.abf', '_Manual_SWDs.csv');

swdManualFile = fullfile(fileDir, [fileBase, '_Manual_SWDs.csv']);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%