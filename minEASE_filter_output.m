function output = minEASE_filter_output (varargin)
%% Filter the combined output from a minEASE output directory
% Usage: output = minEASE_filter_output (varargin)
% Examples:
%       minEASE_filter_output
%       minEASE_filter_output('OutputDir', someDir)
%       minEASE_filter_output('OutputFiles', {file1, file2})
%       minEASE_filter_output('OutputDir', someDir, ...
%                               'OutputFiles', {file1, file2})
% Outputs:
%       output      - an output structure array containing the fields:
%                   siMs
%                   nSamples
%                   prevSweepsDuration
%                   outputLabel
%                   eventInfo
%                   eventClass
%                   fileIdentifier
%                   directionLabel
%                   eventInfoFiltered
%                   eventClassFiltered
%                   isCheckedFiltered
%                   interEventIntervals
%                       
% Arguments:
%       varargin    - 'OutputFiles': output file names from a minEASE run
%                   must be a valid file or a cell array of valid files
%                       under outputDir
%                   default == use the default in minEASE_load_output.m
%                   - 'OutputDir': output directory from a minEASE run
%                   must be a valid directory'
%                   default == use the default in minEASE_load_output.m
%                   - 'ClassesToInclude': classes of events to include
%                   must be a positive integer vector
%                   default == 1:3 (Type I~III PSCs)
%                   - 'TimeWindow': time window to look for events in ms
%                   must be a numeric vector
%                   default == [min(timeVector), max(timeVector)]
%
% Requires:
%       cd/minEASE_filter_events.m
%       cd/minEASE_load_output.m
%       cd/dlmwrite_with_header.m
%
% Used by:
%       cd/ZG_extract_all_IEIs.m

% File History:
% 2018-07-29 Created by Adam Lu
% 2018-08-01 Appended classesToInclude and timeWindow to the output file names
% 2018-08-03 Renamed sweepLabel -> outputLabel
% 2018-08-03 Fixed placement of timeWindowStr in the code
% 2018-08-03 Now uses prevSweepsDuration for timeStart
% 2018-08-03 Added isCheckedFiltered
% 2018-08-12 Implement the case when timeUnits is 'ms'
% 

%% Hard-coded constants
filteredStr = 'ALL_FILTERED';

%% Column assignments for eventInfo:
%   1 = index or time (ms) at event breakpoint
%   2 = index or time (ms) at event peak
%   3 = data value at event breakpoint
%   4 = data value at event peak
%   5 = amplitude of the event
%   6 = 0-100% rise time (samples or ms)
%   7 = 10-90% rise time (samples or ms)
%   8 = inter-event interval from this event peak 
%       to next event peak (samples or ms)
%   9 = interstimulus interval from this event peak 
%       to next event breakpoint (samples or ms)
%   10 = 50% decay time (samples or ms)
%   11 = "full decay" time (samples or ms):
%           time to return within noiseLevel 
%           of breakpoint value
outputColumnHeaderSamples = {'Breakpoint Time (samples)', ...
                            'Peak Time (samples)', ...
                            'Breakpoint Value (pA)', ...
                            'Peak Value (pA)', ...
                            'Peak Amplitude (pA)', ...
                            '0-100% Rise Time (samples)', ...
                            '10-90% Rise Time (samples)', ...
                            'Peak to Peak Interval (samples)', ...
                            'Peak to Breakpoint Interval (samples)', ...
                            '50% Decay Time (samples)', ...
                            'Full Decay Time (samples)', ...
                            'Event Class', ...
                            'Whether Examined'};
outputColumnHeaderMs = {'Breakpoint Time (ms)', ...
                        'Peak Time (ms)', ...
                        'Breakpoint Value (pA)', ...
                        'Peak Value (pA)', ...
                        'Peak Amplitude (pA)', ...
                        '0-100% Rise Time (ms)', ...
                        '10-90% Rise Time (ms)', ...
                        'Peak to Peak Interval (ms)', ...
                        'Peak to Breakpoint Interval (ms)', ...
                        '50% Decay Time (ms)', ...
                        'Full Decay Time (ms)', ...
                        'Event Class', ...
                        'Whether Examined'};

%% Default values for optional arguments
outputFileNamesDefault = {};            % to be set later
outputDirDefault = '';                  % to be set later
classesToIncludeDefault = 1:3;
timeWindowDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputFiles', outputFileNamesDefault, ...
    @(x) assert(ischar(x) || iscell(x) && (all(cellfun(@ischar, x)) ...
        || all(cellfun(@isstring, x))) || isstring(x) , ...
        ['OutputFileNames must be either a string/character array ', ...
            'or a cell array of strings/character arrays!']));
addParameter(iP, 'OutputDir', outputDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'ClassesToInclude', classesToIncludeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'vector'}));
addParameter(iP, 'TimeWindow', timeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, varargin{:});
outputFileNames = iP.Results.OutputFiles;
outputDir = iP.Results.OutputDir;
classesToInclude = iP.Results.ClassesToInclude;
timeWindowUser = iP.Results.TimeWindow;

%% Load the combined output matfile in the directory
% Load everything
%   Note: Default outputDir and outputFileNames will be defined here
output = minEASE_load_output('OutputFiles', outputFileNames, ...
                             'OutputDir', outputDir);

% Count the number of outputs
nFiles = numel(output);

%% Filter the outputs
% Run through all files
for iFile = 1:nFiles
    % Extract information from output structure
    siMs = output(iFile).siMs;
    nSamples = output(iFile).nSamples;
    prevSweepsDuration = output(iFile).prevSweepsDuration;
    outputLabel = output(iFile).outputLabel;
    eventInfo = output(iFile).eventInfo;
    eventClass = output(iFile).eventClass;
    isChecked = output(iFile).isChecked;
    fileIdentifier = output(iFile).fileIdentifier;
    directionLabel = output(iFile).directionLabel;
    timeUnits = output(iFile).timeUnits;

    if isnan(nSamples)
        error('Something isn''t right!');
    end
    
    %% Preparation
    % Get the starting time in milliseconds
    timeStart = prevSweepsDuration;

    % Create a default time window if not provided
    if isempty(timeWindowUser)
        timeWindow = [timeStart, timeStart + nSamples * siMs];
    else
        timeWindow = timeWindowUser;
    end

    % Create strings for outputs
    eventClassStr = ['eventClass', sscanf(num2str(classesToInclude), '%s')];
    timeWindowStr = ['timeWindow', strjoin(strsplit(num2str(timeWindow)), 'to')];

    % Create filtered output file names
    outputFileBase = [fileIdentifier, '_', directionLabel, '_', ...
                        filteredStr, '_', eventClassStr, '_', ...
                        timeWindowStr, '_', timeUnits];
    filteredOutputMatFileName = [outputFileBase, '.mat'];
    filteredOutputCsvFileName = [outputFileBase, '.csv'];
    filteredOutputCsvWHeaderFileName = [outputFileBase, '_w_header.csv'];
    matFileName = fullfile(outputDir, filteredOutputMatFileName);
    csvFileName = fullfile(outputDir, filteredOutputCsvFileName);
    csvWHeaderFileName = fullfile(outputDir, filteredOutputCsvWHeaderFileName);

    %% Filter the events
    % Filter events
    [eventInfoFiltered, eventClassFiltered, ...
        isCheckedFiltered, interEventIntervals] = ...
        minEASE_filter_events(eventInfo, eventClass, isChecked, ...
                        nSamples, siMs, timeStart, ...
                        'ClassesToInclude', classesToInclude, ...
                        'TimeWindow', timeWindow, ...
                        'EventInfoTimeUnits', timeUnits);

    %% Save and output results
    % Save filtered results as matfile
    save(matFileName, 'eventInfoFiltered', 'eventClassFiltered', ...
                    'isCheckedFiltered', 'interEventIntervals', ...
                    'siMs', 'nSamples', 'prevSweepsDuration', ...
                    'outputLabel', '-v7.3');

    % Save filtered results as a comma-separated value file
    outputMatrix = [eventInfoFiltered, eventClassFiltered, isCheckedFiltered];
    csvwrite(csvFileName, outputMatrix);

    % Save filtered results as a comma-separated value file with header
    if strcmp(timeUnits, 'samples')
        dlmwrite_with_header(csvWHeaderFileName, outputMatrix, ...
                                'ColumnHeader', outputColumnHeaderSamples);                   
    elseif strcmp(timeUnits, 'ms')
        dlmwrite_with_header(csvWHeaderFileName, outputMatrix, ...
                                'ColumnHeader', outputColumnHeaderMs);                   
    else
        error('Not implemented yet!');
    end

    % Place in output structure
    output(iFile).eventInfoFiltered = eventInfoFiltered;
    output(iFile).eventClassFiltered = eventClassFiltered;
    output(iFile).isCheckedFiltered = isCheckedFiltered;
    output(iFile).interEventIntervals = interEventIntervals;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

save(fullfile(outputDir, filteredOutputMatFileName), '-struct', 'output', '-v7.3');

% Create filtered output matfile name
filteredOutputMatFileName = [fileIdentifier, '_', directionLabel, '_', ...
                             filteredOutputSuffix, '.mat'];
filteredOutputCsvFileName = [fileIdentifier, '_', directionLabel, '_', ...
                             filteredOutputSuffix, '.csv'];

timeStart = prevSweepsDuration - siMs;

outputDirDefault = pwd;

output = minEASE_load_output('OutputFiles', outputFileNames, ...
                             'OutputDir', outputDir, 'LoadMode', 'combined');

if ~isfolder(outputDir)
    fprintf('%s does not exist or is not readable!\n', outputDir);
    return;
end

timeUnits = 'samples';

%}
