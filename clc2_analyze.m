%% Analyzes all CLC2 data
%
% Requires: 
%       cd/parse_all_multiunit.m

% File History:
% 2019-08-06 Created by Adam Lu
% 

%% Hard-coded parameters
parentDir = '/media/adamX/CLC2/data/blinded';
dirsToAnalyze = {'drug-clean', 'control-clean'};

% For compute_default_signal2noise.m
relSnrThres2Max = 0.1;

% For detect_spikes_multiunit.m
filtFreq = [100, 1000];
minDelayMs = 25;

% For compute_spike_density.m
binWidthMs = 10;                % use a bin width of 10 ms by default
resolutionMs = 5;

% For compute_spike_histogram.m
minBurstLengthMs = 20;          % bursts must be at least 20 ms by default
maxInterBurstIntervalMs = 1000; % bursts are no more than 
                                %   1 second apart by default
minSpikeRateInBurstHz = 100;    % bursts must have a spike rate of 
                                %   at least 100 Hz by default

% For compute_autocorrelogram.m
filterWidthMs = 100;
minRelProm = 0.02;

% For compute_phase_average.m
nSweepsLastOfPhase = 10;         % select from last 10 values by default
nSweepsToAverage = 5;            % select 5 values by default
maxRange2Mean = 40;              % range is not more than 40% of mean 
                                 %   by default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run through all directories
for iDir = 1:numel(dirsToAnalyze)
    % Get the current directory to analyze
    dirThis = fullfile(parentDir, dirsToAnalyze{i});

    % Parse all slices in this directory
    parse_all_multiunit('Directory', dirThis, ...
                        'SaveMatFlag', true, ...
                        'PlotSpikeDetection', true, ...
                        'PlotMeasures', true, ...
                        'PlotSpikeDensity', true, ...
                        'RelSnrThres2Max', relSnrThres2Max, ...
                        'FiltFreq', filtFreq, ...
                        'MinDelayMs', minDelayMs, ...
                        'BinWidthMs', binWidthMs, ...
                        'ResolutionMs', resolutionMs, ...
                        'MinBurstLengthMs', minBurstLengthMs, ...
                        'MaxInterBurstIntervalMs', maxInterBurstIntervalMs, ...
                        'MinSpikeRateInBurstHz', minSpikeRateInBurstHz, ...
                        'FilterWidthMs', filterWidthMs, ...
                        'MinRelProm', minRelProm, ...
                        'NSweepsLastOfPhase', nSweepsLastOfPhase, ...
                        'NSweepsToAverage', nSweepsToAverage, ...
                        'MaxRange2Mean', maxRange2Mean);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%