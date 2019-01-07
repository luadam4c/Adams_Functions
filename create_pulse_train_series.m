function waveform = create_pulse_train_series (varargin)
%% Creates a pulse train series (a theta burst stimulation by default)
% Usage: waveform = create_pulse_train_series (varargin)
% Explanation:
%       Creates a pulse train series (defaults to a theta burst stimulation)
%           i.e., a specific number of sweeps (with a specific duration) 
%                   each containing:
%                   a series (with a specific duration aftera delay) of 
%                   trains (with a specific frequency and duration) of 
%                   pulses (with a specific amplitude, frequency and duration)
%
% Side Effects:
%       (1) Plots and saves the waveform
%       (2) Writes the waveform to a spreadsheet file (default .dat)
%
% Example(s):
%       waveform = create_pulse_train_series;
%       waveform = create_pulse_train_series('PulseAmplitude', 1);
%       waveform = create_pulse_train_series('guiFlag', false, ...
%                                            'plotFlag', false, ...
%                                            'saveFlag', false);
%
% Outputs:
%       waveform    - a waveform for the pulse train series
%                   specified as a column vector
%
% Arguments:    
%       varargin    - 'GuiFlag': whether to open a dialog box to check params
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotFlag': whether to plot the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveFlag': whether to save the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'OutFolder': directory to output files and figures
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SheetType': sheet type; 
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'dat'
%                   - 'FigTypes': figure type(s) for saving
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - 'SamplingRate': sampling rate in Hz
%                   must be a positive scalar
%                   default == 10000 Hz
%                   - 'PulseAmplitude': pulse amplitude
%                   must be a numeric scalar
%                   default == 5
%                   - 'PulseFrequency': pulse frequency in Hz
%                   must be a positive scalar
%                   default == 100 Hz
%                   - 'TrainFrequency': train frequency in Hz
%                   must be a positive scalar
%                   default == 5 Hz
%                   - 'PulseDuration': pulse duration in ms
%                   must be a positive scalar
%                   default == 1 ms
%                   - 'TrainDuration': train duration in ms
%                   must be a positive scalar
%                   default == 40 ms
%                   - 'SeriesDelay': series delay in ms
%                   must be a positive scalar
%                   default == 1000 ms
%                   - 'SeriesDuration': series duration in ms
%                   must be a positive scalar
%                   default == 2000 ms
%                   - 'TotalDuration': total duration in ms
%                   must be a positive scalar
%                   default == 10000 ms
%                   - 'NSweeps': total number of sweeps
%                   must be a positive integer scalar
%                   default == 10
%                   
%
% Requires:
%       cd/create_pulse.m
%       cd/create_waveform_train.m
%       cd/isfigtype.m
%       cd/issheettype.m
%       cd/print_or_show_message.m
%       cd/save_all_figtypes.m

% File History:
% 2018-08-08 Created by Adam Lu
% 2018-08-09 Added seriesDelay and NSweeps
% 2018-08-09 Now defaults the pulse amplitude to 5
% 2018-08-10 Added an input dialog box
% 2018-12-13 Pulled code out to create_pulse.m

% Hard-coded constants
MS_PER_S = 1000;                % 1000 ms per second

% Hard-coded parameters
dialogDimensions = [1, 60];

%% Default values for optional arguments
guiFlagDefault = true;          % opens a dialog box to check params by default
plotFlagDefault = true;         % plot the pulse train series by default
saveFlagDefault = true;         % save the pulse train series by default
outFolderDefault = '';          % default directory to output spreadsheet file
sheetTypeDefault = 'dat';       % default spreadsheet type
figTypesDefault = 'png';        % default figure type(s) for saving
samplingRateDefault = 10000;    % default sampling rate is 10 kHz
pulseAmplitudeDefault = 5;      % default pulse amplitude is 5 V 
                                %   (industry standard for an ON signal)
pulseFrequencyDefault = 100;    % default pulse frequency is 100 Hz
trainFrequencyDefault = 5;      % default train frequency is 5 Hz
pulseDurationDefault = 1;       % default pulse duration is 4 ms
trainDurationDefault = 40;      % default train duration is 40 ms
seriesDelayDefault = 1000;      % default series delay is 1 second
seriesDurationDefault = 2000;   % default series duration is 2 seconds
totalDurationDefault = 10000;   % default total duration is 10 seconds
nSweepsDefault = 10;            % default sweep number is 10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'GuiFlag', guiFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'SamplingRate', samplingRateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'PulseAmplitude', pulseAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PulseFrequency', pulseFrequencyDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TrainFrequency', trainFrequencyDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'PulseDuration', pulseDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TrainDuration', trainDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'SeriesDelay', seriesDelayDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'SeriesDuration', seriesDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TotalDuration', totalDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'NSweeps', nSweepsDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, varargin{:});
guiFlag = iP.Results.GuiFlag;
plotFlag = iP.Results.PlotFlag;
saveFlag = iP.Results.SaveFlag;
outFolder = iP.Results.OutFolder;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
samplingRate = iP.Results.SamplingRate;
pulseAmplitude = iP.Results.PulseAmplitude;
pulseFrequency = iP.Results.PulseFrequency;
trainFrequency = iP.Results.TrainFrequency;
pulseDuration = iP.Results.PulseDuration;
trainDuration = iP.Results.TrainDuration;
seriesDelay = iP.Results.SeriesDelay;
seriesDuration = iP.Results.SeriesDuration;
totalDuration = iP.Results.TotalDuration;
nSweeps = iP.Results.NSweeps;

%% Create a parameters dialog box for the user to modify/check
% Prepare for parameters dialog box
prompt = {'Sampling rate in Hz:', ...
            'Pulse amplitude:', ...
            'Pulse frequency in Hz:', ...
            'Train frequency in Hz:', ...
            'Pulse duration in ms:', ...
            'Train duration in ms:', ...
            'Series delay in ms:', ...
            'Series duration in ms:', ...
            'Total duration in ms:', ...
            'Total number of sweeps:'};
dialogTitle = 'Parameters for the pulse train series';
dimensions = repmat(dialogDimensions, numel(prompt), 1);
defaultParams = {num2str(samplingRate), ...
                num2str(pulseAmplitude), ...
                num2str(pulseFrequency), ...
                num2str(trainFrequency), ...
                num2str(pulseDuration), ...
                num2str(trainDuration), ...
                num2str(seriesDelay), ...
                num2str(seriesDuration), ...
                num2str(totalDuration), ...
                num2str(nSweeps)};

% Open parameters dialog box:
if guiFlag
    % Inputs need to be checked
    inputsValid = false;
else
    % Inputs already valid
    inputsValid = true;
end
mTitle = 'Preference invalid';  % for a message box if needed
icon = 'help';
while ~inputsValid
    % Open input dialog box
    params = inputdlg(prompt, dialogTitle, dimensions, defaultParams, 'on');

    % If the user closes it, exit the program
    if isempty(params)
        print_or_show_message('Action cancelled ...', ...
                              'MTitle', 'Create Pulse Train Series Warning', ...
                              'Icon', 'warn');
        return;
    end

    % Read in user inputs
    samplingRate = str2double(params{1});
    pulseAmplitude = str2double(params{2});
    pulseFrequency = str2double(params{3});
    trainFrequency = str2double(params{4});
    pulseDuration = str2double(params{5});
    trainDuration = str2double(params{6});
    seriesDelay = str2double(params{7});
    seriesDuration = str2double(params{8});
    totalDuration = str2double(params{9});
    nSweeps = str2double(params{10});

    % Update defaults shown in the new dialog box if unsuccessful
    defaultParams = params;

    if isnan(samplingRate) || samplingRate <= 0
        msg = 'Sampling Rate must be a positive scalar!';
    elseif isnan(pulseAmplitude)
        msg = 'Pulse Amplitude must be a numeric scalar!';
    elseif isnan(pulseFrequency) || pulseFrequency <= 0
        msg = 'Pulse Frequency must be a positive scalar!';
    elseif isnan(trainFrequency) || trainFrequency <= 0
        msg = 'Train Frequency must be a positive scalar!';
    elseif isnan(pulseDuration) || pulseDuration <= 0
        msg = 'Pulse Duration must be a positive scalar!';
    elseif isnan(trainDuration) || trainDuration <= 0
        msg = 'Train Duration must be a positive scalar!';
    elseif isnan(seriesDelay) || seriesDelay <= 0
        msg = 'Series Delay must be a positive scalar!';
    elseif isnan(seriesDuration) || seriesDuration <= 0
        msg = 'Series Duration must be a positive scalar!';
    elseif isnan(totalDuration) || totalDuration <= 0
        msg = 'Total Duration must be a positive scalar!';
    elseif isnan(nSweeps) || nSweeps <= 0 || round(nSweeps) ~= nSweeps
        msg = 'NSweeps must be a positive integer!';
    elseif pulseDuration > floor(MS_PER_S / pulseFrequency)
        msg = 'Pulse Duration too large or Pulse Frequency too large!';
    elseif pulseDuration > trainDuration
        msg = 'Pulse Duration cannot be greater than Train Duration!';
    elseif trainDuration > floor(MS_PER_S / trainFrequency)
        msg = 'Train Duration too large or Train Frequency too large!';
    elseif trainDuration > seriesDuration
        msg = 'Train Duration cannot be greater than Series Duration!';
    elseif seriesDelay + seriesDuration > totalDuration
        msg = ['The sum of Series Delay and Series Duration ', ...
                'cannot be greater than Total Duration!'];
    else                        % all inputs are valid
        msg = '';
    end

    % Check whether there was an invalid input
    if isempty(msg)
        % Exit while loop
        inputsValid = true;
    else
        % Show error message and wait for user to close it
        %   before reloading session preferences dialog box
        uiwait(msgbox(msg, mTitle, 'modal'));    
    end
end

%% Preparation
% Compute the sampling interval in ms
siMs = (1 / samplingRate) * MS_PER_S;

% Create a file base
fileBase = ['pulseTrainSeries', '_', num2str(pulseAmplitude), 'V', ...
            num2str(pulseDuration), 'ms', num2str(pulseFrequency), 'Hz', ...
            'Pulse', '_', num2str(trainDuration), 'ms', ...
            num2str(trainFrequency), 'Hz', 'Train', '_', ...
            num2str(seriesDuration), 'ms', 'Series', '_', ...
            num2str(seriesDelay), 'ms', 'Delay', '_', ...
            num2str(totalDuration), 'ms', 'Swp', '_', ...
            num2str(samplingRate/1000), 'kHz'];


%% Create a pulse train
% Create a pulse
pulse = create_pulse('SamplingInterval', siMs, 'PulseDelay', 0, ...
                    'PulseDuration', pulseDuration, ...
                    'PulseAmplitude', pulseAmplitude);

% Create a pulse train from the pulse
pulseTrain = create_waveform_train(pulse, pulseFrequency, trainDuration, ...
                                        'SamplingRate', samplingRate, ...
                                        'plotFlag', false, 'saveFlag', false);
if isempty(pulseTrain)
    waveform = [];
    return;
end

%% Create a pulse train series from the pulse train
trainSeries = create_waveform_train (pulseTrain, trainFrequency, ...
                                seriesDuration, 'SamplingRate', samplingRate, ...
                                'plotFlag', false, 'saveFlag', false);
if isempty(trainSeries)
    waveform = [];
    return;
end
seriesLength = length(trainSeries);

%%  Create the waveform
% Convert delay and total duration to samples
seriesDelaySamples = floor(seriesDelay / siMs);
totalDurationSamples = floor(totalDuration / siMs);

% Create the waveform
waveform = zeros(totalDurationSamples, 1);
seriesIndices = seriesDelaySamples + (1:seriesLength);
waveform(seriesIndices) = trainSeries;

% Create a time vector
timeVec = (1:totalDurationSamples)' * siMs;

% Place the time and waveform vector in the same array
%   with nSweeps repetition counts
waveformArray = [timeVec, repmat(waveform, 1, nSweeps)];

%% Plot the waveform
if plotFlag
    % Plot the waveform
    h = figure;
    clf(h)
    plot(timeVec, waveform);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title(['Stimulus waveform for ', fileBase], 'Interpreter', 'none');

    % Create the full path to the fileBase
    fileBasePath = fullfile(outFolder, fileBase);

    % Save the figure
    save_all_figtypes(h, fileBasePath, figTypes);
end

%% Write output
if saveFlag
    % Append the number of sweeps to the file name
    %   if greater than 1
    if nSweeps > 1
        sheetBase = [fileBase, '_', num2str(nSweeps), 'Swps'];
    else
        sheetBase = fileBase;
    end

    % Create the full path to the spreadsheet file
    sheetPath = fullfile(outFolder, [sheetBase, '.', sheetType]);

    % Convert the waveform to a table
    waveformTable = array2table(waveformArray);

    % Write the waveform to a spreadsheet file
    writetable(waveformTable, sheetPath, ...
                'WriteVariableNames', false, ...
                'Delimiter', '\t');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute the pulse train duration
trainDurationSamples = length(pulseTrain);

% Compute the repeat interval in ms
repeatIntervalMs = (1 / repeatRate) * MS_PER_S;

% Compute the repeat interval in samples
repeatIntervalSamples = floor(repeatIntervalMs / siMs);

% Compute the number of repeats
nRepeats = floor(repeatDurationSamples / repeatIntervalSamples);

% Create a single repeat
singleRepeat = zeros(repeatIntervalSamples, 1);
singleRepeat(1:trainDurationSamples) = pulseTrain;

% Create all repeats
allRepeats = repmat(singleRepeat, nRepeats, 1);

% Convert pulse duration to samples
pulseDurationSamples = floor(pulseDuration / siMs);

% Create a pulse
pulse = ones(pulseDurationSamples, 1) * pulseAmplitude;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%