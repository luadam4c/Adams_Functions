function [eventInfo, eventClass, isChecked, siMs] = ...
                minEASE_combine_events (varargin)
%% Convert info from Excel file into a structure for GUI
% Usage: [eventInfo, eventClass, isChecked, siMs] = ...
%               minEASE_combine_events (varargin)
% Outputs:
%   TODO
%
% Arguments:
%   TODO
%       varargin    - 'Folder': folder containing files to combine
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'ExpLabel': experiment label for files to combine
%                   must be a string scalar or a character vector
%                   default == '' 
%                   - 'OutputLabel': experiment label for output file names
%                   must be a string scalar or a character vector
%                   default == expLabel if provided and folder name otherwise
%                   - 'TimeUnits': units used for time
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'ms'        - milliseconds
%                       'samples'   - samples
%                   default == 'samples'
%                   - 'MessageMode': how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%                   - 'Verbose': whether to print to standard output
%                                   regardless of message mode
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
% Requires:
%       cd/minEASE_extract_from_output_filename.m
%       /home/Matlab/Adams_Functions/dlmwrite_with_header.m
%       /home/Matlab/Adams_Functions/find_in_strings.m
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%
% Used by:
%       cd/minEASE.m
%
% File History:
% 2017-07-24 Created by AL
% 2017-07-25 AL - Added 'units'; now saves two versions -- 
%                   one in samples and one in ms
% 2017-07-25 AL - Now saves 'expLabel' & 'units' too
% 2017-08-01 AL - Now returns if there is nothing to combine
% 2017-10-15 AL - Added success message
% 2017-10-16 AL - Now allows eventInfo to be empty
% 2018-01-28 AL - Added isdeployed
% 2018-02-02 AL - Swp*.mat -> Swp*_output.mat
% 2018-02-02 AL - Added showMessage as an optional parameter-value pair argument
% 2018-02-02 AL - Now uses print_or_show_message.m for output
% 2018-02-07 MD - Changed usage of print_or_show_message()
% 2018-02-14 AL - Added 'output' to output files
% 2018-02-27 AL - Changed showMessages to messageMode with possible values:
%                   'wait', 'show', 'none'
% 2018-03-02 MD - Defined verbose parameter for print_or_show_message
% 2018-05-29 AL - Fixed usage of dlmwrite_with_header()
% 2018-08-02 AL - Made 'TimeUnits' an optional parameter
% 2018-08-02 AL - Made 'ExpLabel' an optional parameter and 
%                   changed it to use '' as the default
% 2018-08-02 AL - Made 'OutputLabel' an optional parameter and 
%                   changed it to use the provided ExpLabel or 
%                   the output directory name as the default
% 2018-08-02 AL - Made 'Folder' an optional parameter
% 2018-08-02 AL - Changed eventInfo -> eventInfoThis; allEventInfo -> eventInfo
%                   eventClass -> eventClassThis; allEventClass -> eventClass
% 2018-08-02 AL - Now stores nSamples, siMs, prevSweepsDuration, sweepLabel
% 2018-08-03 AL - Updated sweepLabelToMatch to include underscores on each side
% TODO: Improve performance by storing and using nEvents for each 
%           individual matfile
% TODO: Check if siMs are the same for each cell
% TODO: Check whether outputs already exist, and deal with them
% TODO: Check whether the number of sweeps available equals the number of files
%

%% Hard-coded parameters
validMessageModes = {'wait', 'show', 'none'};
validTimeUnits = {'samples', 'ms'};
individualOutputSuffix = '_output.mat';

%% Constants to be consistent with find_directional_events.m
IDXBREAK_COLNUM      = 1;
IDXPEAK_COLNUM       = 2;
VALBREAK_COLNUM      = 3;
VALPEAK_COLNUM       = 4;
EVENTAMP_COLNUM      = 5;
TOTALRISE_COLNUM     = 6;
TENNINETYRISE_COLNUM = 7;
IEI_COLNUM           = 8;
ISI_COLNUM           = 9;
HALFDECAY_COLNUM     = 10;
FULLDECAY_COLNUM     = 11;

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
folderDefault = pwd;            % use the present working directory by default
expLabelDefault = '';           % disregard any experiment label by default
outputLabelDefault = '';        % (will be changed later)
timeUnitsDefault = 'samples';   % use samples for the time units by default
messageModeDefault = 'none';    % print to standard output by default
verboseDefault = false;         % default: Program does not print message
                                %   even if message box is shown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Folder', folderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ExpLabel', expLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutputLabel', outputLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
folder = iP.Results.Folder;
expLabel = iP.Results.ExpLabel;
outputLabel = iP.Results.OutputLabel;
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);
verbose = iP.Results.Verbose;


%% Preparation
tic;

% Decide on output cell array header
switch timeUnits
case 'samples'
    outputColumnHeader = outputColumnHeaderSamples;
case 'ms'
    outputColumnHeader = outputColumnHeaderMs;
otherwise
    error('Problem with code!\n');
end

% TODO: Check whether outputs already exist, and deal with them

% Find all matfiles in this folder with this expLabel (may be empty)
files = dir(fullfile(folder, [expLabel, '*Swp*', individualOutputSuffix]));
fileNames = {files.name};

% Count the number of output matfiles (number of sweeps)
nFiles = numel(fileNames);

% Return if there is no event information to combine
if nFiles == 0
    % Show warning message
    message = 'There is no output file to combine!';
    mTitle = 'Combine Event Info Warning';
    icon = 'warn';
    print_or_show_message(message, 'MessageMode', messageMode, ...
                            'MTitle', mTitle, 'Icon', icon, 'Verbose', verbose);
    return;
end

% Get the directionLabel from the first output file name
infoStruct = minEASE_extract_from_output_filename(fileNames{1});
directionLabel = infoStruct.directionLabel;

% Decide on the output label
[~, folderName, ~] = fileparts(folder);
if isempty(outputLabel)
    if ~isempty(expLabel)
        outputLabel = expLabel;
    else
        outputLabel = [folderName, '_', directionLabel];
    end
end

% Decide on the output file names
matFileName = ...
    fullfile(folder, sprintf('%s_ALL_output_%s.mat', outputLabel, timeUnits));
csvFileName = ...
    fullfile(folder, sprintf('%s_ALL_output_%s.csv', outputLabel, timeUnits));
csvFileNameWHeader = ...
    fullfile(folder, sprintf('%s_ALL_output_%s_w_header.csv', ...
                            outputLabel, timeUnits));

% Initialize outputs
eventInfo = [];
eventClass = [];
isChecked = [];
siMs = [];
prevSweepsSamples = 0;          % previous sweep durations 
                                %   for this experiment (samples)
prevSweepsDuration = 0;         % previous sweep durations 
                                %   for this experiment (ms)
prevSweepsDurationInitial = 0;

toc;

%% Combine eventInfo etc.
tic;

% Assume the number of sweeps available equals the number of files
%   TODO: Check this first
nSweeps = nFiles;

for iSweep = 1:nSweeps
    % Get the current sweep label to match
    sweepLabelToMatch = ['_Swp', num2str(iSweep), '_'];

    % Find the file with the file name containing the sweep label
    [indFile, matchedFile] = ...
        find_in_strings(sweepLabelToMatch, fileNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true); 

    % Make sure such a file exists
    if isempty(indFile)
        % Show warning message
        message = {sprintf(['Missing output file containing the ', ...
                            'substrings ''%s'' and ''%s''! '], ...
                            expLabel, sweepLabelToMatch), ...
                    'Cannot combine eventInfo of sweeps!'};
        mTitle = 'Combine Event Info Warning';
        icon = 'warn';
        print_or_show_message(message, 'MessageMode', messageMode, ...
                                'MTitle', mTitle, 'Icon', icon, ...
                                'Verbose', verbose);
        return;
    end

    % Make sure there is no more than one such file
    if length(indFile) > 1
        % Show warning message
        message = {sprintf(['There is more than one file containing the ', ...
                            'substrings ''%s'' and ''%s''! '], ...
                            expLabel, sweepLabelToMatch), ...
                    'Cannot combine eventInfo of sweeps!'};
        mTitle = 'Combine Event Info Warning';
        icon = 'warn';
        print_or_show_message(message, 'MessageMode', messageMode, ...
                                'MTitle', mTitle, 'Icon', icon, ...
                                'Verbose', verbose);
        return;
    end

    % Make sure the direction label is the same as the first file
    infoStruct = minEASE_extract_from_output_filename(matchedFile);
    directionLabelThis = infoStruct.directionLabel;
    if ~strcmp(directionLabel, directionLabelThis)
        % Show warning message
        message = {sprintf(['The direction label for this file is %s, ', ...
                            'but %s is expected!'], ...
                            directionLabelThis, directionLabel), ...
                    'Cannot combine eventInfo of sweeps!'};
        mTitle = 'Combine Event Info Warning';
        icon = 'warn';
        print_or_show_message(message, 'MessageMode', messageMode, ...
                                'MTitle', mTitle, 'Icon', icon, ...
                                'Verbose', verbose);
        return;
    end

    % Extract eventInfo, eventClass, isChecked from mat file
    m = matfile(fullfile(folder, matchedFile));
    eventInfoThis = m.eventInfo;
    eventClassThis = m.eventClass;
    isCheckedThis = m.isChecked;
    siMs = m.siMs;                  % sampling interval (ms)
    nSamplesThis = m.nSamples;          % number of samples

    % If this is the first file, change the previous sweeps duration
    %   to whatever is stored
    if iSweep == 1 && isfield(m, 'prevSweepsDuration')
        prevSweepsDuration = m.prevSweepsDuration;
        prevSweepsDurationInitial = prevSweepsDuration;
    end


    % Compute current sweep duration
    sweepDuration = siMs * nSamplesThis;    % sweep duration in ms

    if ~isempty(eventInfoThis)
        if strcmpi(timeUnits, 'samples')
            % Convert all 'indices' in eventInfoThis to index over the entire experiment
            eventInfoThis(:, IDXBREAK_COLNUM:IDXPEAK_COLNUM) = ...
                prevSweepsSamples + eventInfoThis(:, IDXBREAK_COLNUM:IDXPEAK_COLNUM);
        elseif strcmpi(timeUnits, 'ms')
            % Convert all 'indices' in eventInfoThis to actual time (ms) 
            %   over the entire experiment
            eventInfoThis(:, IDXBREAK_COLNUM:IDXPEAK_COLNUM) = ...
                prevSweepsDuration + ...
                    eventInfoThis(:, IDXBREAK_COLNUM:IDXPEAK_COLNUM) * siMs;

            % Convert all 'times' in eventInfoThis to ms
            eventInfoThis(:, TOTALRISE_COLNUM:FULLDECAY_COLNUM) = ...
                eventInfoThis(:, TOTALRISE_COLNUM:FULLDECAY_COLNUM) * siMs;
        end
    end
 
    % Add to eventInfo, eventClass, isChecked
    eventInfo = [eventInfo; eventInfoThis];
    eventClass = [eventClass; eventClassThis];
    isChecked = [isChecked; isCheckedThis];

    % Update duration of previous sweeps (ms)
    prevSweepsSamples = prevSweepsSamples + nSamplesThis;
    prevSweepsDuration = prevSweepsDuration + sweepDuration;
end

% Record the starting time of the first sweep
prevSweepsDuration = prevSweepsDurationInitial;

% Record the total number of samples
nSamples = prevSweepsSamples;

toc;

%% Save outputs
tic;

% Save eventInfo, eventClass, isChecked as a matfile
save(matFileName, 'eventInfo', 'eventClass', 'isChecked', 'timeUnits', ...
                    'nSamples', 'siMs', 'prevSweepsDuration', ...
                    'expLabel', 'outputLabel', '-v7.3');

% Save eventInfo, eventClass, isChecked as a csv file
%   Note: xlswrite does not work on fishfish
%   Note: csvwrite can only write numeric arrays
outputMatrix = [eventInfo, eventClass, isChecked];
csvwrite(csvFileName, outputMatrix);

% Save eventInfo, eventClass, isChecked with header as a csv file
dlmwrite_with_header(csvFileNameWHeader, outputMatrix, ...
                        'ColumnHeader', outputColumnHeader);

%% Show success message
message = {sprintf('Eventinfo successfully combined using time in %s!!', ...
                    timeUnits), ...
            'You will love it! ^^'};
mTitle = 'Combine Event Info Success';
icon = 'none';
% TODO: implement icon
print_or_show_message(message, 'MessageMode', messageMode, ...
                        'MTitle', mTitle, 'Icon', icon, 'Verbose', verbose);

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

fprintf('There is no output file to combine!');
fprintf('Missing %s! Cannot combine eventInfo of sweeps!\n\n', ...
            matfilename);
fprintf('Eventinfo successfully combined!!\n\n')

%                   - 'ShowMessage': whether to show messages in messages boxes
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
showMessageDefault  = false;        % print to standard output by default
addParameter(iP, 'ShowMessage', showMessageDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
showMessage = iP.Results.ShowMessage;
        if showMessage
            print_or_show_message(message, 'MessageMode', 'show', ...
                                    'MTitle', mTitle, 'Icon', icon);
        else
            print_or_show_message(message, 'MessageMode', 'none', ...
                                    'MTitle', mTitle, 'Icon', icon);
        end

files = dir(fullfile(folder, [expLabel, '_Swp*output.mat']));

matFileName = ...
    fullfile(folder, sprintf('%s_ALL_output_%s.mat', expLabel, timeUnits));
csvFileName = ...
    fullfile(folder, sprintf('%s_ALL_output_%s.csv', expLabel, timeUnits));
csvFileNameWHeader = ...
    fullfile(folder, sprintf('%s_ALL_output_%s_w_header.csv', expLabel, timeUnits));

matfilename = fullfile(folder, ...
                [expLabel, '_Swp', num2str(iSweep), '_output.mat']);
m = matfile(matfilename);
if isempty(whos(m))
    % Show warning message
    message = sprintf('Missing %s! Cannot combine eventInfo of sweeps!', ...
                matfilename);

sweepLabel = ['Swp', num2str(iSweep)];

for iSweep = 76:nSweeps
sweepLabelToMatch = ['_Swp', num2str(files.name{iSweep})), '_'];

%}
