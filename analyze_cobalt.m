%% Clear workspace
clear all force hidden
close all force hidden

%% Hard coded parameters
lowpassNpoles = 8;                 	% order of the ButterWorth filter

channelsOfInterest = [13, 15, 17, 19];
% lowpassCutoffs = [0.5, 4, 8, 14, 30];
% highpassCutoffs = [3, 7, 13, 30, 80];
lowpassCutoffs = [0.5, 4, 8, 14, 30] * 100;
highpassCutoffs = [3, 7, 13, 30, 80] * 100;
bandName = {'delta', 'theta', 'alpha', 'beta', 'gamma'};

%% Construct paths
scriptsPath = pwd;
[homeDir, ~, ~] = fileparts(scriptsPath);
[parentDir, ~, ~] = fileparts(homeDir);
eeglabDir = fullfile(parentDir, 'eeglab14_1_2b');
dataDir = fullfile(homeDir, 'data');
outputDir = fullfile(homeDir, 'output');

%% Add functions to search path
addpath(eeglabDir);

%% Preparation
nChannels = length(channelsOfInterest);

%% EXTRACT DATA
% Change the directory to the data directory
cd(dataDir);

% Read in RHD2000 data
%   These variables created will be used:
%       amplitier_data, t_amplifier, filename
read_Intan_RHD2000_file;

% Change the directory back to the scripts directory
cd(scriptsPath);

% Check if data exists
if exist('amplifier_data', 'var') ~= 1
    error('amplifier_data does not exist!');
end

% Get the base of the file name
% fileBase = strrep(filename, '.rhd', '');
[~, fileBase, ~] = fileparts(filename);

% Read in the time vector in seconds
timeVec = t_amplifier';

% Compute the sampling interval in seconds
si = timeVec(2) - timeVec(1);

% Compute the number of samples
nSamples = length(timeVec);

% Transpose the amplifier data so that each column is a signal
amplifierDataFlipped = amplifier_data';

% Extract data from channels of interest
data = amplifierDataFlipped(:, channelsOfInterest);

% Save data
dataFileName = fullfile(outputDir, [fileBase, '.mat']);
save(dataFileName, 'data', '-v7.3');

%% Exploratory Analysis
% Perform analysis with Signal Analyzer
% signalAnalyzer(data)

%% Bandpass filter the signal
% Compute the number of bands to filter
nBands = length(lowpassCutoffs);

% Construct file suffices
fileSuffices = cell(nBands, 1);
figTitles = cell(nBands, 1);
parfor iBand = 1:nBands
    fileSuffices{iBand} = ['_', bandName{iBand}, '_band'];
    figTitles{iBand} = ['Data filtered ', ...
                        num2str(lowpassCutoffs(iBand)), ...
                        '~', num2str(highpassCutoffs(iBand)), ...
                        ' Hz (', bandName{iBand}, ' band)'];
end

% Preallocate bandpass-filtered data in a cell array
dataFiltered = cell(nBands, 1);

% Filter the data with a zero-phase bandpass Butterworth filter
%   with 3 dB cutoff frequencies [lowpassCutoff, highpassCutOff]
%   and order lowpassNpoles each way (or 2*lowpassNpoles in total)
for iBand = 1:nBands
    % Get the current bandpass cutoffs
    lowpassCutoff = lowpassCutoffs(iBand);
    highpassCutoff = highpassCutoffs(iBand);
    
    % Find the normalized cutoff frequency Wn = fc/(fs/2), 
    %   where fs = sampling frequency (Hz) = 1/si 
    %   and fs/2 is the Nyquist frequency
    Wn = [lowpassCutoff, highpassCutoff] * 2 * si;
                            % normalized cutoff frequency (half-cycles/sample)

    % Find the transfer function coefficients of a lowpass Butterworth filter
    %   with order npoles and normalized cutoff frequency Wn
    [numeratorCoeff, denominatorCoeff] = butter(lowpassNpoles, Wn, 'bandpass');

    % Check the order of the filter
    orderFilter = filtord(numeratorCoeff, denominatorCoeff);
    if nSamples <= 3 * orderFilter
        error(['Not enough data points to apply a ', ...
                'Butterworth filter of order %d twice!\n'], ...
                orderFilter);
    end

    % Bandpass-filter data twice (forward & reverse directions)
    dataFiltered{iBand} = filtfilt(numeratorCoeff, denominatorCoeff, data);
end

%% Plot figures

% Create a figure for plotting raw data
% Create a file name
figName = fullfile(outputDir, [fileBase, '_raw_data', '.jpg']);
h = figure(10000);
clf(h)
plot_signals(timeVec, data, 'Title', 'Raw data');
saveas(h, figName);
%close(h)

for iBand = 1:nBands
    figName = fullfile(outputDir, [fileBase, fileSuffices{iBand}, '.jpg']);
    h = figure(10000 + iBand);
    clf(h)
    plot_signals(timeVec, dataFiltered{iBand}, 'Title', figTitles{iBand});
    saveas(h, figName);
end

%%
%{
OLD CODE

exampleFile = fullfile(dataDir, 'SmartboxRecording_20180620-112301.rhd');

% Change the directory to the output directory
cd(outputDir);

%}