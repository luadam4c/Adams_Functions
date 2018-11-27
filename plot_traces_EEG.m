function plot_traces_EEG (abfFileName, varargin)
%% Plots EEG traces from a .abf file
% Usage: plot_traces_EEG (abfFileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       plot_traces_EEG('WAGS04_30_2018_cage1.abf');
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
%                   - 'OutFolder': the name of the directory in which 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == a subdirectory named by {fileName}_traces in pwd
%                   - 'AtfPath': path to .atf file
%                   must be a TODO
%                   default == TODO
%                   - 'AssystPath': path to Assyst.txt file
%                   must be a TODO
%                   default == TODO
%                   - 'SayliPath': path to Sayli_SWDs.csv file
%                   must be a TODO
%                   default == TODO
%                   
% Requires:
%       cd/all_files.m
%       cd/check_dir.m
%       cd/parse_abf.m
%       cd/parse_assyst_swd.m
%       cd/parse_atf_swd.m
%       cd/plot_traces_abf.m
%
% Used by:
%       cd/plot_EEG.m

% File History:
%   2016-09-XX - Created
%   2018-11-21 - Renamed spikewavedetection -> plot_EEG
%   2018-11-21 - Updated to use plot_traces_abf.m
%   2018-11-26 - Removed first argument and detect in current directory
% 

%% Default values for optional arguments
verboseDefault = true;
outFolderDefault = '';      % set later
atfPathDefault = '';        % set later
assystPathDefault = '';     % set later
sayliPathDefault = '';      % set later

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
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'AtfPath', atfPathDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'AssystPath', assystPathDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SayliPath', sayliPathDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, abfFileName, varargin{:});
verbose = iP.Results.Verbose;
outFolder = iP.Results.OutFolder;
atfPath = iP.Results.AtfPath;
assystPath = iP.Results.AssystPath;
sayliPath = iP.Results.SayliPath;

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

% Decide on the output folder
if isempty(outFolder)
    outFolder = fullfile(fileDir, strcat(fileBase, '_traces'));
end
if verbose
    fprintf('Outfolder is %s ...\n', outFolder);
end

% Check if needed output directory exist
check_dir(outFolder, 'Verbose', verbose);

% Look for .atf files if not provided
if isempty(atfPath)
    % Try to look for a .atf file containing 
    %   the file base in the same directory
    [~, atfPaths] = all_files('Directory', fileDir, 'Keyword', fileBase, ...
                                'Extension', '.atf');

    % If not successful, try to look for a .atf file containing 
    %   the file base in the 'atffiles' directory
    if isempty(atfPaths)
        [~, atfPaths] = ...
            all_files('Directory', fullfile(fileDir, 'atffiles'), ...
                        'Keyword', fileBase, 'Extension', '.atf');
    end

    % If there is more than one file, one choose the first one
    if isempty(atfPaths)
        atfPath = '';
    elseif iscell(atfPaths)
        atfPath = atfPaths{1};
    else
        atfPath = atfPaths;
    end
end

% Look for Assyst .txt files if not provided
if isempty(assystPath)
    % Try to look for a Assyst.txt file containing 
    %   the file base in the same directory
    [~, assystPaths] = all_files('Directory', fileDir, 'Keyword', fileBase, ...
                                'Suffix', 'Assyst', 'Extension', '.txt');

    % If not successful, try to look for a Assyst.txt file containing 
    %   the file base in the 'Assyst' directory
    if isempty(assystPaths)
        [~, assystPaths] = ...
            all_files('Directory', fullfile(fileDir, 'Assyst'), ...
                        'Keyword', fileBase, ...
                        'Suffix', 'Assyst', 'Extension', '.txt');
    end

    % If there is more than one file, one choose the first one
    if isempty(assystPaths)
        assystPath = '';
    elseif iscell(assystPaths)
        assystPath = assystPaths{1};
    else
        assystPath = assystPaths;
    end
end

% Look for Sayli SWD files if not provided
if isempty(sayliPath)
    % Try to look for a Sayli_SWDs.csv file containing 
    %   the file base in the same directory
    [~, sayliPaths] = all_files('Directory', fileDir, 'Keyword', fileBase, ...
                                'Suffix', 'Sayli_SWDs', 'Extension', '.csv');

    % If not successful, try to look for a Sayli_SWDs.csv file containing 
    %   the file base in the 'Sayli' directory
    if isempty(sayliPaths)
        [~, sayliPaths] = ...
            all_files('Directory', fullfile(fileDir, 'Sayli'), ...
                        'Keyword', fileBase, ...
                        'Suffix', 'Sayli_SWDs', 'Extension', '.csv');
    end

    % If there is more than one file, one choose the first one
    if isempty(sayliPaths)
        sayliPath = '';
    elseif iscell(sayliPaths)
        sayliPath = sayliPaths{1};
    else
        sayliPath = sayliPaths;
    end
end

%% Do the job
% Plot full traces if not already done so
plot_traces_abf(abfPath, 'Verbose', verbose, 'OverWrite', false, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
                'OutFolder', outFolder);

% Plot manual SWD regions if an atf file exists
if isfile(atfPath)
    % Parse the SWDs from the .atf file
    swdManualTable = parse_atf_swd(atfPath);

    % Get the recorded .abf path
    abfPathRecordedCell = swdManualTable{1, 'abfPath'};
    abfPathRecorded = abfPathRecordedCell{1};

    % Check if the abfPaths are the same as the provided abfPath
    if ~strcmp(abfPath, abfPathRecorded)
        fprintf('The recorded path %s does not match %s!!\n', ...
                abfPathRecorded, abfPath);
    else
        % TODO: Pull out to function plot_swds.m
        % Create an output folder for manual SWD traces
        outFolderSWD = fullfile(outFolder, 'Manual_SWDs');

        % Plot each detected SWD separately
        plot_traces_abf(abfPath, 'Verbose', verbose, 'OverWrite', false, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
                'TimeStart', swdManualTable.startTime, ...
                'TimeEnd', swdManualTable.endTime, ...
                'OutFolder', outFolderSWD);
    end
end

% Plot Assyst SWD regions if an Assyst.txt file exists
if isfile(assystPath)
    % Parse the SWDs from the Assyst.txt file
    swdAssystTable = parse_assyst_swd(assystPath);

    % Get the recorded .abf path
    abfPathRecordedCell = swdAssystTable{1, 'abfPath'};
    abfPathRecorded = abfPathRecordedCell{1};

    % Check if the abfPaths are the same as the provided abfPath
    if ~strcmp(abfPath, abfPathRecorded)
        fprintf('The recorded path %s does not match %s!!\n', ...
                abfPathRecorded, abfPath);
    else
        % TODO: Pull out to function plot_swds.m
        % Create an output folder for Assyst SWD traces
        outFolderSWD = fullfile(outFolder, 'Assyst_SWDs');

        % Plot each detected SWD separately
        plot_traces_abf(abfPath, 'Verbose', verbose, 'OverWrite', false, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
                'TimeStart', swdAssystTable.startTime, ...
                'TimeEnd', swdAssystTable.endTime, ...
                'OutFolder', outFolderSWD);
    end
end

% Plot Sayli SWD regions if an Sayli_SWDs.csv file exists
if isfile(sayliPath)
    % Parse the SWDs from the Assyst.txt file
    swdSayliTable = readtable(sayliPath);

    % TODO: Pull out to function plot_swds.m
    % Create an output folder for Sayli SWD traces
    outFolderSWD = fullfile(outFolder, 'Sayli_SWDs');

    % Plot each detected SWD separately
    plot_traces_abf(abfPath, 'Verbose', verbose, 'OverWrite', false, ...
            'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
            'TimeStart', swdSayliTable.startTimes, ...
            'TimeEnd', swdSayliTable.finishTimes, ...
            'OutFolder', outFolderSWD);
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
if isempty(atfPath)
    atfPath = replace(abfPath, '.abf', '_Manual_SWDs.csv');
end

% Construct full path to abf file
[abfPath, pathExists] = construct_and_check_abfpath(abfFileName);
if ~pathExists
    return
end

atfPath = replace(abfPath, '.abf', '_Manual_SWDs.csv');

atfPath = fullfile(fileDir, [fileBase, '_Manual_SWDs.csv']);

swdAssystTable = parse_assyst_swd(assystPath, 'OutFolder', outFolder);
swdManualTable = parse_atf_swd(atfPath, 'OutFolder', outFolder);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%