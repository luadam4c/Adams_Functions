function waveformTrain = create_waveform_train (waveform, frequency, totalDuration, varargin)
%% Creates a waveform train from a waveform, a frequency, and a total duration
% Usage: waveformTrain = create_waveform_train (waveform, frequency, totalDuration, varargin)
% Explanation:
%       Creates a waveform train based on a waveform, the waveform frequency
%           in Hz and the total duration in milliseconds
%
% Side Effects:
%       (1) Plots and saves the waveform
%       (2) Writes the waveform to a spreadsheet file (default .dat)
%
% Example(s):
%       waveform = ones(10, 1) * 5;
%       waveformTrain = create_waveform_train (waveform, 100, 40)
%
% Outputs:
%       waveformTrain  - a train of waveforms
%                   specified as a numeric column vector
%
% Arguments:    
%       waveform    - a waveform
%                   must be a numeric vector
%       frequency   - waveform frequency in Hz
%                   must be a positive scalar
%       totalDuration   - total train duration in ms
%                   must be a positive scalar
%       varargin    - 'SamplingRate': sampling rate in Hz
%                   must be a positive scalar
%                   default == 10000 Hz
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
%                   - 'PlotFlag': whether to plot the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveFlag': whether to save the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       /home/Matlab/Adams_Functions/isfigtype.m
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%       /home/Matlab/Adams_Functions/save_all_figtypes.m
%
% Used by:    
%       /home/Matlab/Adams_Functions/create_pulse_train_series.m

% File History:
% 2018-08-08 Created by Adam Lu


% Hard-coded constants
MS_PER_S = 1000;                % 1000 ms per second

% Hard-coded parameters

%% Default values for optional arguments
samplingRateDefault = 10000;    % default sampling rate is 10 kHz
outFolderDefault = '';          % default directory to output spreadsheet file
sheetTypeDefault = 'dat';       % default spreadsheet type
figTypesDefault = 'png';        % default figure type(s) for saving
plotFlagDefault = true;         % plot the pulse train series by default
saveFlagDefault = true;         % save the pulse train series by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'waveform', ...             % a waveform
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'frequency', ...            % waveform frequency in Hz
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addRequired(iP, 'totalDuration', ...        % total train duration in ms
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SamplingRate', samplingRateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
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
parse(iP, waveform, frequency, totalDuration, varargin{:});
samplingRate = iP.Results.SamplingRate;
outFolder = iP.Results.OutFolder;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
plotFlag = iP.Results.PlotFlag;
saveFlag = iP.Results.SaveFlag;

%% Preparation
% Make sure the waveform is a column
waveform = waveform(:);

% Compute the sampling interval in ms
siMs = (1 / samplingRate) * MS_PER_S;

% Create a file base
fileBase = ['waveformTrain', '_', num2str(totalDuration), 'ms', ...
            num2str(frequency), 'Hz', 'Train'];

% Count the waveform duration in samples
waveformDurationSamples = length(waveform);

% Count the total number of samples in the train
totalDurationSamples = totalDuration / siMs;

% Compute the waveform interval in ms
waveformIntervalMs = (1 / frequency) * MS_PER_S;

% Compute the waveform interval in samples
waveformIntervalSamples = floor(waveformIntervalMs / siMs);

% Return with error message if the waveform is longer than the interval
if waveformDurationSamples > waveformIntervalSamples
    waveformTrain = [];
    message = {['The waveform duration cannot be longer ', ...
                    'than the waveform interval.'], ...
                'Please decrease the waveform frequency!'};
    mTitle = 'Create Waveform Error';
    icon = 'error';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                        'MessageMode', 'show', 'Verbose', true, ...
                        'CreateMode', 'non-modal');
    return
end

%% Create the waveformTrain
% Compute the number of waveforms
nWaveforms = floor(totalDurationSamples / waveformIntervalSamples);

% Create a single waveform
singleWaveform = zeros(waveformIntervalSamples, 1);
singleWaveform(1:waveformDurationSamples) = waveform;

% Create all waveforms
allWaveforms = repmat(singleWaveform, nWaveforms, 1);

% Find the length of all waveforms
allPulsesLength = length(allWaveforms);

% Create the waveform train
waveformTrain = zeros(totalDurationSamples, 1);
waveformTrain(1:allPulsesLength) = allWaveforms;

% Create a time vector
timeVec = (1:totalDurationSamples)' * siMs;

% Place the time and waveform vector in the same array
waveformArray = [timeVec, waveformTrain];

% Plot the waveform
if plotFlag
    % Plot the waveform
    h = figure;
    clf(h)
    plot(timeVec, waveformTrain);
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
    % Create the full path to the spreadsheet file
    sheetPath = fullfile(outFolder, [fileBase, '.', sheetType]);

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

% Convert the waveform to a table
waveformTable = array2table(waveformTrain);

%}