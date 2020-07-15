function [swdManualTable, swdManualCsvFile] = ...
                parse_atf_swd (originalEventFile, varargin)
%% Parse spike-wave-discharge (SWD) event info from .atf file or converted .csv file
% Usage: [swdManualTable, swdManualCsvFile] = ...
%               parse_atf_swd (originalEventFile, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       swdManualTable = parse_atf_swd('WAGS04_30_2018_cage3_Manual_SWDs.atf');
%
% Outputs:
%       swdManualTable      - a table of spike-wave discharge event info
%                           specified as a 2D table
%
% Arguments:
%       originalEventFile   - original event file, could be .atf or 
%                               converted .csv
%                           must be a string scalar or a character vector
%       varargin    - 'TraceFileName': Name of the corresponding trace file(s)
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - 'OutFolder': directory to output swd table file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == same as location of originalEventFile
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   - 'ParseFeatures': whether to parse SWD features
%                                       such as peakFrequency
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/argfun.m
%       cd/array_fun.m
%       cd/atf2sheet.m
%       cd/construct_and_check_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/create_label_from_sequence.m
%       cd/create_time_vectors.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/force_column_cell.m
%       cd/is_overlapping.m
%       cd/issheettype.m
%       cd/match_dimensions.m
%       cd/parse_abf.m
%       cd/read_data_atf.m
%       cd/read_lines_from_file.m
%       cd/sscanf_full.m
%       cd/union_over_cells.m
%
% Used by:
%       cd/parse_all_swds.m
%       cd/parse_psd.m
%       cd/plot_traces_EEG.m

% File History:
% 2018-11-21 Created by Adam Lu
% 2018-12-26 Added 'SheetType' as an optional argument
% 2019-09-08 Now uses the trace file name as the basis 
%               for constructing sheet file name
% 2019-09-09 Updated the construction of trace file paths
% 2019-09-24 Added check for overlapping windows
% 2020-06-26 Added 'ParseFeatures' as an optional argument
% 2020-07-15 Now generates an output file even there is no event recorded

%% Hard-coded constants
MS_PER_S = 1000;
N_LINES_TO_SKIP = 2;            % scored atf files have two irrelevant lines

%% Hard-coded parameters
varNames = {'startTime', 'endTime', 'duration', 'tracePath', 'pathExists'};

%% Default values for optional arguments
traceFileNameDefault = '';      % set later
outFolderDefault = '';          % set later
sheetTypeDefault = 'csv';       % default spreadsheet type
parseFeaturesDefault = false;

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
addRequired(iP, 'originalEventFile', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TraceFileName', traceFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'ParseFeatures', parseFeaturesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, originalEventFile, varargin{:});
traceFileName = iP.Results.TraceFileName;
outFolder = iP.Results.OutFolder;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
parseFeatures = iP.Results.ParseFeatures;

%% Preparation
% Decide what file type the first input is
if regexpi(originalEventFile, '.atf$')
    atfFile = originalEventFile;
    atfCsvFile = '';
elseif regexpi(originalEventFile, '_atf.csv$')
    atfFile = '';
    atfCsvFile = originalEventFile;
else
    atfFile = '';
    atfCsvFile = '';
end

% Get the file directory
fileDir = fileparts(originalEventFile);

% Decide on the output folder
if isempty(outFolder)
    outFolder = fileDir;
end

% Compute the number of header lines in the scored ATF file
nHeaderLinesAtf = N_LINES_TO_SKIP + 1;

%% Do the job
% Initialize outputs
swdManualTable = table.empty;
swdManualCsvFile = '';

% Read the table from the file
if ~isfile(atfFile) && ~isfile(atfCsvFile)
    % Do nothing
    return
elseif isfile(atfCsvFile)
    % Display warning if atf file also provided
    if isfile(atfFile)
        fprintf(['Table with be read from the csv file ', ...
                'instead of the atf file!\n']);
    end

    % Read in the SWD manual table from the converted csv file
    atfTable = readtable(atfCsvFile);
elseif isfile(atfFile)
    % Make sure it is not an ATF trace file
    isScoredAtf = is_scored_atf(atfFile);
    if ~isScoredAtf
        fprintf('%s is not a scored ATF file!\n', atfFile);
        return;
    end

    % Read in the SWD manual table and print to a csv file
    [atfTable, atfCsvFile] = atf2sheet(atfFile, 'SheetType', sheetType);
    fprintf('%s created!\n', atfCsvFile);
end

% Make sure there is an event recorded
if height(atfTable) == 0
    startTime = [];
    endTime = [];
    duration = [];
    tracePath = {};
    pathExists = [];

    % Construct manual SWD table csv file
    swdManualCsvFile = ...
        replace(atfCsvFile, '_scored_atf.csv', ['_Manual_SWDs.', sheetType]);
else
    % Get the first channel name
    firstSignalName = atfTable{1, 'Signal'};

    % Check whether each row is the same as firstSignalName
    isFirstSignal = strcmp(atfTable.Signal, firstSignalName);

    % Compute the number of windows
    nWindows = sum(isFirstSignal);

    % Compute the number of rows for each window
    switch nWindows
         case 0
            nRowsPerWindow = NaN;
         otherwise
            nRowsPerWindow = round(height(atfTable) / nWindows);
    end 

    % Restrict to the entries for the first channel only
    swdManualTableOfInterest = atfTable(isFirstSignal, :);

    % Get the start and end times in ms
    startTimesMs = swdManualTableOfInterest.Time1_ms_;
    endTimesMs = swdManualTableOfInterest.Time2_ms_;

    % Convert to seconds
    startTime = startTimesMs / MS_PER_S;
    endTime = endTimesMs / MS_PER_S;

    % Make sure none of the windows overlap
    [isOverlapping, ~, indOverlapPrev] = ...
        is_overlapping(transpose([startTime, endTime]));
    if isOverlapping
        toSubtract = (nRowsPerWindow - 1):-1:0;
        rowsOverlapPrevCell = ...
            arrayfun(@(x) nHeaderLinesAtf + x * nRowsPerWindow - toSubtract, ...
                    indOverlapPrev, 'UniformOutput', false);
        rowsOverlapPrev = union_over_cells(rowsOverlapPrevCell);

        fprintf(['The file %s cannot be parsed because the following ', ...
                    'window numbers overlap with the next one:\n'], ...
                    originalEventFile);
        fprintf('\t%s\n', create_label_from_sequence(indOverlapPrev));
        fprintf('Please remove these lines from %s:\n', originalEventFile);
        fprintf('\t%s\n', create_label_from_sequence(rowsOverlapPrev));
        fprintf('\n');
        return
    end

    % Compute duration
    duration = endTime - startTime;

    % If not provided, read in the trace file names
    if isempty(traceFileName)
        % Get the .abf file name for each SWD
        traceFileName = swdManualTableOfInterest.FileName;
    end

    % Construct full path to original data file
    %   Note: Must be (or copied to) the same directory as the original event file
    [tracePath, pathExists] = ...
        construct_and_check_fullpath(traceFileName, 'Directory', fileDir);

    % Force as a column cell array
    tracePath = force_column_cell(tracePath);

    % Make sure the dimensions match up
    [tracePath, pathExists] = ...
        argfun(@(x) match_dimensions(x, size(startTime)), tracePath, pathExists);

    % Extract the file base of the trace file
    traceFileBase = extract_fileparts(traceFileName{1}, 'base');

    % Extract the file extension of the trace file
    traceFileExt = extract_fileparts(traceFileName{1}, 'ext');

    % Construct manual SWD table csv file
    swdManualCsvFile = ...
        fullfile(outFolder, [traceFileBase, '_Manual_SWDs.', sheetType]);    
end

%% Parse SWD features if requested
if parseFeatures
    if height(atfTable) == 0
        peakFrequency = [];
        secondFrequency = [];
        thirdFrequency = [];
        swdFrequency = [];
    else
        [peakFrequency, secondFrequency, thirdFrequency, swdFrequency] = ...
            parse_swd_features(startTime, endTime, tracePath, pathExists);
    end
end

%% Correct the start and end times if the data comes from an ATF file
if height(atfTable) ~= 0
    if strcmpi(traceFileExt, '.atf')
        % Note: The scored atf file itself also has the info, 
        %   but at a weird location

        % If trace file doesn't exist, return error
        if ~isfile(tracePath{1})
            error('The source ATF file %s does not exist!', tracePath{1});
        end

        % Read the sweep start time line
        headerLine = read_lines_from_file(tracePath{1}, 'MaxNum', 1, ...
                                    'Keyword', 'SweepStartTimesMS');
        if isempty(headerLine)
            error('headerLine not found in %s', tracePath{1});
        end

        % Read in the start time of the trace
        traceStartTimeMs = sscanf_full(headerLine, '%g');

        % Extract the start time of the trace
        traceStartTimeSeconds = traceStartTimeMs / MS_PER_S;
    else
        traceStartTimeSeconds = 0;
    end

    % Modify the start and end times
    startTime = traceStartTimeSeconds + startTime;
    endTime = traceStartTimeSeconds + endTime;
end

%% Output results
% Create a table for the parsed SWDs
swdManualTable = table(startTime, endTime, duration, ...
                        tracePath, pathExists, 'VariableNames', varNames);

% Add features to the table
if parseFeatures
    swdManualTable = ...
        addvars(swdManualTable, peakFrequency, secondFrequency, ...
                thirdFrequency, swdFrequency, 'After', 'duration');
end

% Write the table to a file
writetable(swdManualTable, swdManualCsvFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isScoredAtf = is_scored_atf (atfFile)
%% Determine whether the ATF file is scored

% Read in the second line of the file
line2 = read_lines_from_file(atfFile, 'LineNumber', 2);

% Read in the number of lines to skip following this line before the header
%   is reached
nLinesToSkip = sscanf(line2, '%d', 1);

% This number should be zero for a scored ATF file
isScoredAtf = nLinesToSkip == 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [peakFrequency, secondFrequency, thirdFrequency, swdFrequency] = ...
                parse_swd_features (startTime, endTime, tracePath, pathExists)
%% Extracts SWD features from a trace path

[peakFrequency, secondFrequency, thirdFrequency, swdFrequency] = ...
    array_fun(@parse_features_helper, ...
                num2cell(startTime), num2cell(endTime), ...
                tracePath, num2cell(pathExists));
% [peakFrequency, secondFrequency, thirdFrequency, swdFrequency] = ...
%     cellfun(@parse_features_helper, ...
%                 num2cell(startTime), num2cell(endTime), ...
%                 tracePath, num2cell(pathExists));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [peakFrequency, secondFrequency, thirdFrequency, ...
            swdFrequency] = ...
                parse_features_helper (startTime, endTime, tracePath, pathExists)
%% Extracts SWD features from one SWD

%% Hard-coded parameters
filterWindow = 2;
swdFrequencyTarget = 7;

% Note: For pleth experiments done in the Guyunet lab
signalName = 'WIC#2';

% Note: For pleth experiments done in the Beenhakker lab
signalNumber = 1;

%% Preparation
% Make sure path exists
if ~pathExists
    peakFrequency = NaN;
    return
end

% Construct relative time window from start of file
%   Note: this is the times recorded by the scored ATF files
timeWindowRel = [startTime; endTime];

%% Do the job
% Extract signal containing SWD data
if contains(tracePath, '.abf')
    % Read from the ABF file
    [abfParams, abfData] = parse_abf(tracePath, 'TimeUnits', 's');

    % Read the sampling interval in seconds
    siSeconds = abfParams.siSeconds;

    % Read time vector in seconds
    timeVec = abfData.tVec;

    % Choose data vector
    vVecs = abfData.vVecs;
    traceData = vVecs(:, signalNumber);
elseif contains(tracePath, '.atf')
    % Read from the ATF file
    [atfData, atfParams] = read_data_atf('FilePaths', tracePath);

    % Read the sampling interval in seconds
    siSeconds = atfParams.siSeconds;

    % Read time vector in seconds
    timeVec = atfData.Time;

    % Choose data vector
    traceData = atfData.(signalName);
elseif contains(tracePath, '.txt')
    % Note: This is for old pleth data scored by Katie
    % Read in the matrix data
    if get_matlab_year >= 2019
        textData = readmatrix(tracePath, 'FileType', 'text');
    else
        textData = csvread(tracePath);
    end

    % Read the sampling interval in seconds
    siSeconds = 0.005;

    % Count the number of samples
    nSamples = size(textData, 1);

    % Read time vector in seconds
    timeVec = create_time_vectors(nSamples, 'TimeUnits', 's', ...
                                    'SamplingIntervalSeconds', siSeconds);

    % Choose data vector
    traceData = textData(:, signalNumber);
else
    error('The trace path %s is unrecognized!', tracePath);
end

% Find the relative time vector
timeVecRel = timeVec - timeVec(1);

% Find the SWD time window endpoints
endPoints = find_window_endpoints(timeWindowRel, timeVecRel);

% Restrict to SWD region
dataSwd = extract_subvectors(traceData, 'EndPoints', endPoints);

% Compute the power spectral density over the SWD region
[~, psdData] = parse_psd(dataSwd, 'SamplingFrequency', 1/siSeconds, ...
                                    'FilterWindow', filterWindow);
freqSelected = psdData.freqSelected;

% Extract the peak frequencies of the SWD
peakFrequency = freqSelected(1);
secondFrequency = freqSelected(2);
thirdFrequency = freqSelected(3);

% Choose the one closest to swdFrequencyTarget as the SWD frequency
[~, swdFrequency] = find_closest(freqSelected, swdFrequencyTarget);

% % Compute the spectrogram for the SWD
% [spectData, freqHz, timeInstantsSeconds] = ...
%     compute_spectrogram (traceData, siSeconds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%