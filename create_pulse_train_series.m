function waveform = create_pulse_train_series (varargin)
%% Create a stimulus waveform (a theta burst stimulation by default)
% Usage: waveform = create_pulse_train_series (varargin)
% Explanation:
%       Creates a pulse train series (defaults to a theta burst stimulation)
%           i.e., a sweep (with a specific duration) that contains 
%                   a series (with a specific duration) of 
%                   trains (with a specific frequency and duration) of 
%                   pulses (with a specific amplitude, frequency and duration)
%
% Side Effects:
%       (1) Plots and saves the waveform
%       (2) Writes the waveform to a spreadsheet file (default .dat)
%
% Example(s):
%       waveform = create_pulse_train_series;
%       waveform = create_pulse_train_series('PulseAmplitude', 5);
%
% Outputs:
%       waveform    - a waveform for the pulse train series
%                   specified as a column vector
%
% Arguments:    
%       varargin    - 'SamplingRate': sampling rate in Hz
%                   must be a positive scalar
%                   default == 10000 Hz
%                   - 'PulseAmplitude': pulse amplitude
%                   must be a numeric scalar
%                   default == 1
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
%                   - 'SeriesDuration': series duration in ms
%                   must be a positive scalar
%                   default == 2000 ms
%                   - 'TotalDuration': total duration in ms
%                   must be a positive scalar
%                   default == 10000 ms
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
%                   
%
% Requires:
%       /home/Matlab/Adams_Functions/create_pulse_train.m
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/Adams_Functions/save_all_figtypes.m

% File History:
% 2018-08-08 Created by Adam Lu

% Hard-coded constants
MS_PER_S = 1000;                % 1000 ms per second

% Hard-coded parameters

%% Default values for optional arguments
samplingRateDefault = 10000;    % default sampling rate is 10 kHz
pulseAmplitudeDefault = 1;      % default pulse amplitude is 1 V
pulseFrequencyDefault = 100;    % default pulse frequency is 100 Hz
trainFrequencyDefault = 5;      % default train frequency is 5 Hz
pulseDurationDefault = 1;       % default pulse duration is 1 ms
trainDurationDefault = 40;      % default train duration is 40 ms
seriesDurationDefault = 2000;   % default series duration is 2 seconds
totalDurationDefault = 10000;   % default total duration is 10 seconds
outFolderDefault = '';          % default directory to output spreadsheet file
sheetTypeDefault = 'dat';       % default spreadsheet type
figTypesDefault = 'png';        % default figure type(s) for saving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
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
addParameter(iP, 'SeriesDuration', seriesDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TotalDuration', totalDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
samplingRate = iP.Results.SamplingRate;
pulseAmplitude = iP.Results.PulseAmplitude;
pulseFrequency = iP.Results.PulseFrequency;
trainFrequency = iP.Results.TrainFrequency;
pulseDuration = iP.Results.PulseDuration;
trainDuration = iP.Results.TrainDuration;
seriesDuration = iP.Results.SeriesDuration;
totalDuration = iP.Results.TotalDuration;
outFolder = iP.Results.OutFolder;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Preparation
% Compute the sampling interval in ms
siMs = (1 / samplingRate) * MS_PER_S;

% Create a file base
fileBase = ['pulseTrainSeries', '_', num2str(pulseAmplitude), 'V', ...
            num2str(pulseDuration), 'ms', num2str(pulseFrequency), 'Hz', ...
            'Pulse', '_', num2str(trainDuration), 'ms', ...
            num2str(trainFrequency), 'Hz', 'Train', '_', ...
            num2str(seriesDuration), 'ms', 'Series', '_', ...
            num2str(totalDuration), 'ms', 'Sweep'];


%% Create a pulse train
% Convert pulse duration to samples
pulseDurationSamples = floor(pulseDuration / siMs);

% Create a pulse
pulse = ones(pulseDurationSamples, 1) * pulseAmplitude;

% Create a pulse train from the pulse
pulseTrain = create_pulse_train(pulse, pulseFrequency, trainDuration, ...
                                        'SamplingRate', samplingRate);

%% Create a pulse train series from the pulse train
trainSeries = create_pulse_train (pulseTrain, trainFrequency, ...
                                seriesDuration, 'SamplingRate', samplingRate);
seriesLength = length(trainSeries);

%%  Create the waveform
% Count the number of samples
totalDurationSamples = floor(totalDuration / siMs);

% Create the waveform
waveform = zeros(totalDurationSamples, 1);
waveform(1:seriesLength) = trainSeries;

%% Plot the waveform
% Create a time vector
timeVec = (1:totalDurationSamples)' * siMs;

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

%% Write output
% Create the full path to the spreadsheet file
sheetPath = fullfile(outFolder, [fileBase, '.', sheetType]);

% Convert the waveform to a table
waveformTable = array2table(waveform);

% Write the waveform to a spreadsheet file
writetable(waveformTable, sheetPath, 'WriteVariableNames', false);

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

%}