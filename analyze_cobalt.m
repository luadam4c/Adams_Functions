%% Clear workspace
clear all force hidden
close all force hidden

%% Flags
dataType = 'mat'; %'rhd';
plotRawFlag = true; %false;
bandpassFilterFlag = true; %false;
plotFilteredFlag = true; %false;
signalAnalyzerFlag = false; %true;
eeglabFlag = false; %true;
crossCorrFlag = true; %false;

%% Hard coded parameters
lowpassNpoles = 4;                   	% order of the ButterWorth filter

channelsOfInterest = [13, 15, 17, 19];
% lowCutoffs = [0.5, 4, 8, 14, 30];
% highCutoffs = [3, 7, 13, 30, 80];
% lowCutoffs = [0.5, 4, 8, 14, 30] * 100;
% highCutoffs = [3, 7, 13, 30, 80] * 100;
% bandName = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
% lowCutoffs = [0, 0.1, 0.2, 0.3, 0.4] * samplingRate;
% highCutoffs = [0.1, 0.2, 0.3, 0.4, 0.5] * samplingRate;
% bandName = {'band1', 'band2', 'band3', 'band4', 'band5'};
lowCutoffs = [0, 300];
highCutoffs = [300, 4000];
bandName = {'0to300', '300to4000'};

%% Construct paths
parentDir = 'C:\Users\Pinn Analysis\Desktop\Shinnosuke';
homeDir = fullfile(parentDir, 'cobalt');
eeglabDir = fullfile(parentDir, 'eeglab14_1_2b');
scriptsPath = fullfile(homeDir, 'scripts');
dataDir = fullfile(homeDir, 'data');
outFolder = fullfile(homeDir, 'output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add functions to search path
addpath(eeglabDir);
addpath(scriptsPath);

%% Preparation
nChannels = length(channelsOfInterest);

%% EXTRACT DATA
% Change the directory to the data directory
cd(dataDir);

switch dataType
    case 'rhd'
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

        % Read in the time vector in seconds
        timeVec = t_amplifier';

        % Compute the sampling interval in seconds
        si = timeVec(2) - timeVec(1);

        % Compute the sampling rate in Hz
        samplingRate = 1/si;
        
        % Compute the number of samples
        nSamples = length(timeVec);

        % Transpose the amplifier data so that each column is a signal
        amplifierDataFlipped = amplifier_data';

        % Extract data from channels of interest
        dataRaw = amplifierDataFlipped(:, channelsOfInterest);

    case 'mat'
        % Open a GUI to select .mat file exported by LabChart
        [file, path, filterindex] = ...
                uigetfile('*.mat', 'Select an MATLAB Data File', ...
                            'MultiSelect', 'off');

        % Exit if unsuccessful
        if (file == 0)
            return;
        end

        % Load file
        filename = fullfile(path, file);
        fileStruct = load(filename);
        
        % Resave any data_block1 as data
        if isfield(fileStruct, 'data_block1')
            dataRaw = fileStruct.data_block1;
        elseif isfield(fileStruct, 'data_block2')
            dataRaw = fileStruct.data_block2;
        else
            fprintf('Cannot reognize any data in this matfile!\n');
            return;
        end

        % Get the folder name
        [folderName, ~, ~] = fileparts(path);

        % Find the sampling rate from the folder name
        numbers = sscanf_full(folderName, '%d');
        mouseNumber = numbers(1);
        samplingRate = numbers(2);
        
        % Transpose so that each channel is a column
        dataRaw = dataRaw';
        
        % Get the sampling interval in seconds
        si = 1/samplingRate;

        % Count the number of samples
        nSamples = size(dataRaw, 1);

        % Create a time vector in seconds for the entire channel
        timeVec = (1:nSamples)' * si;
    otherwise
end

% Count the number of channels
nChannels = size(dataRaw, 2);

% Get the base of the file name
[~, fileBase, ~] = fileparts(filename);

% Save data
dataFileName = fullfile(outFolder, [fileBase, '.mat']);
save(dataFileName, 'dataRaw', '-v7.3');

%% Exploratory Analysis
% Perform analysis with Signal Analyzer
if signalAnalyzerFlag
    signalAnalyzer(dataRaw)
end

% Perform analysis with EEG lab
if eeglabFlag
    dataTransposed = dataRaw';
    eeglab
end
    
%% Bandpass filter the signal
if bandpassFilterFlag
    % Compute the number of bands to filter
    nBands = length(lowCutoffs);

    % Construct file suffices
    fileSuffices = cell(nBands, 1);
    figTitles = cell(nBands, 1);
    parfor iBand = 1:nBands
        fileSuffices{iBand} = ['_', bandName{iBand}, '_band'];
        figTitles{iBand} = ['Data filtered ', ...
                            num2str(lowCutoffs(iBand)), ...
                            '~', num2str(highCutoffs(iBand)), ...
                            ' Hz (', bandName{iBand}, ' band)'];
    end

    % Preallocate bandpass-filtered data in a cell array
    dataFiltered = cell(nBands, 1);

    % Filter the data with a zero-phase bandpass Butterworth filter
    %   with 3 dB cutoff frequencies [lowCutoff, highCutoff]
    %   and order lowpassNpoles each way (or 2*lowpassNpoles in total)
    for iBand = 1:nBands
        % Get the current bandpass cutoffs
        lowCutoff = lowCutoffs(iBand);
        highCutoff = highCutoffs(iBand);

        % Find the normalized cutoff frequency Wn = fc/(fs/2) 
        %   (half-cycles/sample)
        %   where fs = sampling frequency (Hz) = 1/si 
        %   and fs/2 is the Nyquist frequency
        % Then compute the transfer function coefficients 
        %   of a lowpass Butterworth filter
        %   with order npoles and normalized cutoff frequency Wn
        if lowCutoff <= 0
            Wn = highCutoff * 2 * si;
            [numeratorCoeff, denominatorCoeff] = ...
                butter(lowpassNpoles, Wn, 'low');
        elseif highCutoff >= samplingRate / 2
            Wn = lowCutoff * 2 * si;
            [numeratorCoeff, denominatorCoeff] = ...
                butter(lowpassNpoles, Wn, 'high');
        else
            Wn = [lowCutoff, highCutoff] * 2 * si;
            [numeratorCoeff, denominatorCoeff] = ...
                butter(lowpassNpoles, Wn, 'bandpass');
        end

        % Check the order of the filter
        orderFilter = filtord(numeratorCoeff, denominatorCoeff);
        if nSamples <= 3 * orderFilter
            error(['Not enough data points to apply a ', ...
                    'Butterworth filter of order %d twice!\n'], ...
                    orderFilter);
        end

        % Bandpass-filter data twice (forward & reverse directions)
        dataFiltered{iBand} = filtfilt(numeratorCoeff, ...
                                        denominatorCoeff, dataRaw);
    end  
end

%% Cross correlation analysis
if crossCorrFlag
    if bandpassFilterFlag
        for iBand = 1:nBands
            fileBaseCorr = [fileBase, ' band', num2str(lowCutoffs(iBand)), ...
                                      'to', num2str(highCutoffs(iBand))];

            corrProf = ...
                crosscorr_profile(dataFiltered{iBand}, samplingRate, ...
                                    outFolder, fileBaseCorr);
        end
    else
        corrProf = crosscorr_profile(dataRaw, samplingRate, ...
                                    outFolder, fileBase);
    end
    
end

%% Plot figures
if plotRawFlag
    figFolder = fullfile(outFolder, 'raw');
    if exist(figFolder, 'dir') ~= 7
        mkdir(figFolder);
    end
    % Create a figure for plotting raw data
    % Create a file name
    figName = fullfile(figFolder, [fileBase, '_raw_data', '.jpg']);
    h = figure(10000);
    clf(h)
    plot_signals(timeVec, dataRaw, 'Title', 'Raw data');
    saveas(h, figName);
    %close(h)
end

if plotFilteredFlag
    figFolder = fullfile(outFolder, 'bandpassed');
    if exist(figFolder, 'dir') ~= 7
        mkdir(figFolder);
    end
    for iBand = 1:nBands
        figName = fullfile(figFolder, [fileBase, fileSuffices{iBand}, '.jpg']);
        h = figure(10000 + iBand);
        clf(h)
        plot_signals(timeVec, dataFiltered{iBand}, ...
                     'Title', figTitles{iBand});
        saveas(h, figName);
    end
end
    
%%
%{
OLD CODE

exampleFile = fullfile(dataDir, 'SmartboxRecording_20180620-112301.rhd');

% Change the directory to the output directory
cd(outFolder);

[homeDir, ~, ~] = fileparts(scriptsPath);
[parentDir, ~, ~] = fileparts(homeDir);

samplingRate = 10000; % 1000            % sampling rate in Hz


%}