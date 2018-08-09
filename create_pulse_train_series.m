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
%                   - 'plotFlag': whether to plot the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'saveFlag': whether to save the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   
%
% Requires:
%       /home/Matlab/Adams_Functions/create_waveform_train.m
%       /home/Matlab/Adams_Functions/isfigtype.m
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/Adams_Functions/save_all_figtypes.m

% File History:
% 2018-08-08 Created by Adam Lu
% 2018-08-09 Added seriesDelay and NSweeps

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
seriesDelayDefault = 1000;      % default series delay is 1 second
seriesDurationDefault = 2000;   % default series duration is 2 seconds
totalDurationDefault = 10000;   % default total duration is 10 seconds
nSweepsDefault = 1; %10;            % default sweep number is 10
outFolderDefault = '';          % default directory to output spreadsheet file
sheetTypeDefault = 'dat';       % default spreadsheet type
figTypesDefault = 'png';        % default figure type(s) for saving
plotFlagDefault = true;         % plot the pulse train series by default
saveFlagDefault = true;         % save the pulse train series by default

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
addParameter(iP, 'SeriesDelay', seriesDelayDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'SeriesDuration', seriesDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TotalDuration', totalDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'NSweeps', nSweepsDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
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
outFolder = iP.Results.OutFolder;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
plotFlag = iP.Results.PlotFlag;
saveFlag = iP.Results.SaveFlag;

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
            num2str(totalDuration), 'ms', 'Swp'];


%% Create a pulse train
% Convert pulse duration to samples
pulseDurationSamples = floor(pulseDuration / siMs);

% Create a pulse
pulse = ones(pulseDurationSamples, 1) * pulseAmplitude;

% Create a pulse train from the pulse
pulseTrain = create_waveform_train(pulse, pulseFrequency, trainDuration, ...
                                        'SamplingRate', samplingRate, ...
                                        'plotFlag', false, 'saveFlag', false);

%% Create a pulse train series from the pulse train
trainSeries = create_waveform_train (pulseTrain, trainFrequency, ...
                                seriesDuration, 'SamplingRate', samplingRate, ...
                                'plotFlag', false, 'saveFlag', false);
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

%}