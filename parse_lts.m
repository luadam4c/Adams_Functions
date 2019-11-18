function varargout = parse_lts (vVec0s, varargin)
%% Finds, plots and classifies the most likely low-threshold spike (LTS) candidate in a voltage trace
% Usage: [parsedParams, parsedData] = parse_lts (vVec0s, siMs (opt), varargin)
% Explanation:
%       TODO
%
% Examples:
%       sweepName = 'C101210_0006_3';
%       matFilesDir = '/media/adamX/m3ha/data_dclamp/take4/matfiles';
%       [data, sweepInfo] = m3ha_import_raw_traces(sweepName, 'Directory', matFilesDir);
%       vVecs = extract_columns(data, 2);
%       siMs = sweepInfo.siMs;
%       [parsedParams, parsedData] = parse_lts(vVecs, siMs, 'StimStartMs', 1000);
%       [parsedParams, parsedData] = parse_lts(vVecs, siMs, 'StimStartMs', 1000, 'PlotFlag', true);
%       [data2, sweepInfo] = m3ha_import_raw_traces(sweepName, 'Directory', matFilesDir, 'ImportMode', 'active');
%       vVecs2 = extract_columns(data2, 2);
%       siMs = sweepInfo.siMs;
%       [parsedParams, parsedData] = parse_lts(vVecs2, siMs, 'StimStartMs', 1000);
%
% Outputs: 
%       parsedParams    - a table containing the parsed parameters, 
%                           each row corresponding to a vector, with fields:
%                           actVhold
%                           maxNoise
%                           peakTime
%                           peak2ndDer
%                           peakProm
%                           peakWidth
%                           peakClass
%                           spikesPerPeak
%                           ltsPeakTime
%                           ltsPeakValue
%                           maxSlopeTime
%                           maxSlopeValue
%                           burstTime
%                           spikesPerBurst
%                           spikeThreshold
%                           firstSpikeTime
%                           lastSpikeTime
%                           maxSpikeAmp
%                           minSpikeAmp
%                           spikeFrequency
%                           spikeAdaptation
%                           couldHaveMissed
%
% Arguments:
%       vVec0s      - original voltage vector(s) in mV
%                   must be a numeric array
%                       or a cell array of numeric vectors
%       siMs        - (opt) sampling interval in ms
%                           If not provided, 'tVec0s' must be provided
%                   must be a positive vector
%       varargin    - 'Verbose': whether to output parsed results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotFlag': whether to plot traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveMatFlag': whether to save data as mat file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveSheetFlag': whether to save params as a spreadsheet
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'FileBase': base of filename (without extension) 
%                                   corresponding to each vector
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == 'unnamed_1', 'unnamed_2', ...
%                   - 'OutFolder': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Prefix': prefix to prepend to output file names
%                   must be a character array
%                   default == extract_common_prefix(fileBase)
%                   - 'StimStartMs': time of stimulation start (ms)
%                   must be a positive scalar
%                   default == find_first_jump(vVec0s)
%                   - 'MinPeakDelayMs': minimum peak delay (ms)
%                               after the end of the pulse
%                   must be a nonnegative scalar
%                   default == 0 ms
%                   - 'NoiseWindowMsOrMaxNoise': maximum noise in mV 
%                                           or noise window in ms
%                                           if numel == 1, maxNoise;
%                                           if numel == 2, noiseWindowMs; 
%                   must be a numeric vector
%                   default == computed from the range [0, stimStartMs]
%                   - 'SearchWindowMs': window to search for LTS (ms)
%                   must be within range of tVec0
%                   default == [StimStartMs + minPeakDelayMs, 
%                               tVec0(end) - medfiltWindowMs], 
%                                   where medfiltWindowMs is 30 ms
%                   - 'tVec0s': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   default == created from siMs and vVec0s
%                   - 'tVec2s': time vector(s) for resampling
%                   must be a numeric array & within range of tVec0
%                   default == siMsRes*(round(tVec0(1)/siMsRes):round(tVec0(end)/siMsRes))'
%                   - 'vVec1s': median-filtered voltage vector(s)
%                   must be a numeric array with same length as vVec0
%                   default == medianfilter(vVec0, medfiltWindowMs, siMs)
%                   - 'vVec2s': voltage vector(s) after resampling
%                   must be a numeric array with same length as tVec2
%                   default == interp1(tVec0, vVec1, tVec2, 'linear')
%                   - 'vVec3s': moving-average-filtered vVec1s
%                   must be a numeric array with same length as vVec0
%                   default == movingaveragefilter(vVec1, smoothWindowMs, siMs)
%
% Requires:
%       cd/all_file_bases.m
%       cd/argfun.m
%       cd/check_subdir.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/decide_on_filebases.m
%       cd/extract_common_prefix.m
%       cd/extract_elements.m
%       cd/find_first_jump.m
%       cd/find_in_strings.m
%       cd/m3ha_locate_homedir.m
%       cd/match_time_info.m
%       cd/medianfilter.m
%       cd/movingaveragefilter.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-01-13 Adapted from find_LTS.m
% 2019-02-19 Made siMs an optional argument
% 2019-11-17 Added 'SaveMatFlag' as an optional parameter
% 2019-11-17 Added 'SaveSheetFlag' as an optional parameter
% 2019-11-18 Added 'Prefix' as an optional parameter
% 2019-11-18 Changed mean -> nanmean; std -> nanstd

%% Hard-coded parameters
% Subdirectories in outFolder for placing figures
outSubDirs = {'vtraces', 'LTSanalysis', 'burstanalysis', ...
            'vtraces_scaled', 'gray_area_traces', 'LTScouldbemissed'};

% Parameters used for data analysis
medfiltWindowMs = 30;% width in ms for the median filter for spikes (voltage traces)
smoothWindowMs = 30;% width in ms for the moving average filter for finding narrowest voltage peaks
baseWidthMs = 20;   % width in ms for calculating baseline voltage (holding potential)
ltsThr = -0.0023;        % 2nd derivative in V^2/s^2 below which defines an LTS peak
ltsThrAlt = -0.0081823;  % 2nd derivative in V^2/s^2 above which is the "gray area"
spThr = -45;        % Initial amplitude threshold in mV for detecting a spike 
                    % Will be changed later when LTS peak amplitude is found
                    % 2016-10-14 the highest LTS peak without bursts is -43.8867 mV
                    % 2016-10-18 Note: the highest LTS peak without bursts is actually -48.0006 mV
                    %           the previous value was for a spontaneous LTS, but because of median-filtering,
                    %           -45 mV is probably a safer threshold
spThrRelLts = 10;   % Relative amplitude threshold in mV for detecting a spike above an LTS 
                    % 2016-10-19 The smallest relative amplitude of an action potential riding above the LTS
                    %           is probably between 10~11 mV (see B091010_0006_19)
siMsRes = 1;        % resampling interval in ms (1 kHz)
minSp2PkTime = 0;   % minimum time from the first spike to the peak of the LTS (ms)
slopeSpacing = 1;   % spacing in ms used to calculate slope
mafw3Dv = 3;        % voltage change in mV corresponding to the moving average filter window for finding slopes
%mafw3 = 5;         % width in ms for the moving average filter for finding slopes
slopeSegYHalf = 5;  % how much voltage difference (mV) to plot maxslope line segment below maxslope point

% Parameters used for output files
ltsSheetSuffix = '_ltsParams';
ltsMatSuffix = '_ltsData';
sheetType = 'csv';

%% Default values for optional arguments
siMsDefault = [];               % set later
verboseDefault = true;          % print to standard output by default
plotFlagDefault = false;
saveMatFlagDefault = false;     % don't save parsed data by default
saveSheetFlagDefault = true;    % save parsed params by default
fileBaseDefault = {};           % set later
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
prefixDefault = '';             % set later
stimStartMsDefault = [];        % set later
minPeakDelayMsDefault = 0;
noiseWindowMsORmaxNoiseDefault = [];    % set later
searchWindowMsDefault = [];       % set later
tVec0sDefault = [];             % set later
tVec2sDefault = [];             % set later
vVec1sDefault = [];             % set later
vVec2sDefault = [];             % set later
vVec3sDefault = [];             % set later

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
addRequired(iP, 'vVec0s', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec0s must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'siMs', siMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveSheetFlag', saveSheetFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileBase must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'NoiseWindowMsOrMaxNoise', noiseWindowMsORmaxNoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'SearchWindowMs', searchWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'tVec0s', tVec0sDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec0s must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'tVec2s', tVec2sDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec2s must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVec1s', vVec1sDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec1s must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVec2s', vVec2sDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec2s must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVec3s', vVec3sDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec3s must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vVec0s, varargin{:});
siMs = iP.Results.siMs;
verbose = iP.Results.Verbose;
plotFlag = iP.Results.PlotFlag;
saveMatFlag = iP.Results.SaveMatFlag;
saveSheetFlag = iP.Results.SaveSheetFlag;
fileBase = iP.Results.FileBase;
outFolder = iP.Results.OutFolder;
prefix = iP.Results.Prefix;
stimStartMs = iP.Results.StimStartMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
noiseWindowMsORmaxNoise = iP.Results.NoiseWindowMsOrMaxNoise;
searchWindowMs = iP.Results.SearchWindowMs;
tVec0s = iP.Results.tVec0s;
tVec2s = iP.Results.tVec2s;
vVec1s = iP.Results.vVec1s;
vVec2s = iP.Results.vVec2s;
vVec3s = iP.Results.vVec3s;

%% Preparation
% Count the number of vectors
nVectors = count_vectors(vVec0s);

% Count the number of samples for each vector
nSamples = count_samples(vVec0s);

% Create file bases if not provided
fileBase = decide_on_filebases(fileBase, nVectors);

% Decide on prefix if not provided
if isempty(prefix)
    prefix = extract_common_prefix(fileBase);
end

% Display message
if verbose
    if ~isempty(prefix)
        fprintf('ANALYZING voltage traces for %s ...\n', prefix);
    else
        fprintf('ANALYZING voltage traces ...\n');
    end
end

% Find files to override
[fileBasesToOverride, idxMissedLtsByOrder, idxMissedLtsByShape, ...
    idxSpikesPerBurstIncorrect, idxLooksLikeMissedLts, ...
    idxLooksLikeLtsNotByProminence, ...
    idxLooksLikeLtsNotByNarrowness, ...
    idxNoiseInTrace, idxSpontLtsOrBurst, ...
    idxWideLtsCouldBeNoise] = m3ha_find_files_to_override;

% Compute sampling interval(s) and create time vector(s)
if isempty(siMs) && isempty(tVec0s)
    error('One of siMs and tVec0s must be provided!');
else
    [tVec0s, siMs, nSamples] = ...
        match_time_info(tVec0s, siMs, nSamples, 'TimeUnits', 'ms');
end

% Compute the base of the time vector(s)
tBase = extract_elements(tVec0s, 'first') - siMs;

% Detect stimulation start time if not provided
% TODO: Not tested
if isempty(stimStartMs)
    % Force as a cell array
    vVec0sCell = force_column_cell(vVec0s);
    
    % Look for the first deviant
    [~, idxStimStart] = ...
        cellfun(@(x) find_first_jump(x, 'Signal2noise', 3), vVec0sCell);

    % Compute the stimulation start in ms
    stimStartMs = idxStimStart .* siMs;
end

% Decide whether to detect holding potential
if stimStartMs < tBase + baseWidthMs
    fprintf(['stimStartMs must be at least %g ms after tBase ', ...
                'for baseline voltage computation!\n'], baseWidthMs);
    fprintf('Actual holding potential will be NaN!\n');
    computeActVholdFlag = false;
else
    computeActVholdFlag = true;
end

% Decide on whether to compute maximum noise
if numel(noiseWindowMsORmaxNoise) == 1 || ...
        numel(noiseWindowMsORmaxNoise) == nVectors
    maxNoise = noiseWindowMsORmaxNoise;
    noiseWindowMs = NaN(2, nVectors);
    computeMaxNoiseFlag = false;
else
    noiseWindowMs = noiseWindowMsORmaxNoise;
    if isempty(noiseWindowMs)
        noiseWindowMs = [0; stimStartMs];
    end
    maxNoise = NaN(nVectors, 1);
    computeMaxNoiseFlag = true;
end

% Compute the minimum peak time in ms
minPeakTimeMs = stimStartMs + minPeakDelayMs;

% Force vectors to be a column cell array
[tVec0s, tVec2s, vVec0s, vVec1s, vVec2s, vVec3s, ...
    noiseWindowMs, searchWindowMs, fileBase] = ...
    argfun(@force_column_cell, ...
            tVec0s, tVec2s, vVec0s, vVec1s, vVec2s, vVec3s, ...
            noiseWindowMs, searchWindowMs, fileBase);

% Make sure all parameters are column vectors with the same number of elements
[siMs, maxNoise, stimStartMs, minPeakTimeMs, ...
    tVec0s, tVec2s, vVec0s, vVec1s, vVec2s, vVec3s, ...
    noiseWindowMs, searchWindowMs, fileBase] = ...
    argfun(@(x) match_dimensions(x, [nVectors, 1]), ...
            siMs, maxNoise, stimStartMs, minPeakTimeMs, ...
            tVec0s, tVec2s, vVec0s, vVec1s, vVec2s, vVec3s, ...
            noiseWindowMs, searchWindowMs, fileBase);

% Parse all of them in a parfor loop
parsedParamsCell = cell(nVectors, 1);
parsedDataCell = cell(nVectors, 1);
parfor iVec = 1:nVectors
%for iVec = 1:nVectors
    [parsedParamsCell{iVec}, parsedDataCell{iVec}] = ...
        parse_lts_helper(verbose, plotFlag, computeActVholdFlag, ...
            computeMaxNoiseFlag, outFolder, ...
            tVec0s{iVec}, tVec2s{iVec}, vVec0s{iVec}, ...
            vVec1s{iVec}, vVec2s{iVec}, vVec3s{iVec}, ...
            noiseWindowMs{iVec}, searchWindowMs{iVec}, fileBase{iVec}, ...
            siMs(iVec), maxNoise(iVec), ...
            stimStartMs(iVec), minPeakTimeMs(iVec), ...
            fileBasesToOverride, idxMissedLtsByOrder, idxMissedLtsByShape, ...
            idxSpikesPerBurstIncorrect, idxLooksLikeMissedLts, ...
            idxLooksLikeLtsNotByProminence, idxLooksLikeLtsNotByNarrowness, ...
            idxNoiseInTrace, idxSpontLtsOrBurst, idxWideLtsCouldBeNoise, ...
            medfiltWindowMs, smoothWindowMs, baseWidthMs, ltsThr, ltsThrAlt, ...
            spThr, spThrRelLts, siMsRes, minSp2PkTime, slopeSpacing, ...
            mafw3Dv, slopeSegYHalf, outSubDirs);
end

% Convert to a struct array
%   Note: This removes all entries that are empty
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@(x) [x{:}], parsedParamsCell, parsedDataCell);

% Convert to a table
[parsedParams, parsedData] = ...
    argfun(@(x) struct2table(x, 'AsArray', true), ...
            parsedParamsStruct, parsedDataStruct);

%% Save outputs
if saveMatFlag
    % Create a path for the LTS data .mat file
    ltsMatPath = fullfile(outFolder, [prefix, ltsMatSuffix, '.mat']);

    % Save outputs in the .mat file
    save(ltsMatPath, 'parsedParams', 'parsedData', '-v7.3');
end
if saveSheetFlag
    % Create a path for the LTS info spreadsheet file
    ltsSheetPath = fullfile(outFolder, [prefix, ltsSheetSuffix, ...
                                        '.', sheetType]);

    % Save table to the spreadsheet file
    writetable(parsedParams, ltsSheetPath);
end

%% Outputs
varargout{1} = parsedParams;
varargout{2} = parsedData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
        parse_lts_helper(verbose, plotFlag, computeActVholdFlag, ...
            computeMaxNoiseFlag, outFolder, ...
            tVec0, tVec2, vVec0, vVec1, vVec2, vVec3, ...
            noiseWindowMs, searchWindowMs, fileBase, ...
            siMs, maxNoise, stimStartMs, minPeakTimeMs, ...
            fileBasesToOverride, idxMissedLtsByOrder, idxMissedLtsByShape, ...
            idxSpikesPerBurstIncorrect, idxLooksLikeMissedLts, ...
            idxLooksLikeLtsNotByProminence, idxLooksLikeLtsNotByNarrowness, ...
            idxNoiseInTrace, idxSpontLtsOrBurst, idxWideLtsCouldBeNoise, ...
            medfiltWindowMs, smoothWindowMs, baseWidthMs, ltsThr, ltsThrAlt, ...
            spThr, spThrRelLts, siMsRes, minSp2PkTime, slopeSpacing, ...
            mafw3Dv, slopeSegYHalf, outSubDirs)


%% Preparation for each trace
% Check relationships between arguments
if ~isequal(numel(tVec0), numel(vVec0))
    error('Time and voltage vectors do not have the same length!');
end

% Create default search window
if isempty(searchWindowMs)
    searchWindowMs = [minPeakTimeMs, tVec0(end) - medfiltWindowMs];    
end

% Display message
if verbose
    fprintf('FINDING LTSs for %s ...\n', fileBase);
    fprintf('Sampling interval == %g ms\n', siMs);
end

%% Reorganize data if needed
% Median filter voltage traces to get rid of spikes
if isempty(vVec1)
    vVec1 = medianfilter(vVec0, medfiltWindowMs, siMs);
end

% Resample all traces at 1 kHz if needed
if isempty(tVec2)
    tVec2 = siMsRes * (round(tVec0(1)/siMsRes):round(tVec0(end)/siMsRes))';
end

% Interpolate voltage vector if needed
if isempty(vVec2)
    vVec2 = interp1(tVec0, vVec1, tVec2, 'linear');
end

% Moving-average-filter median-filtered traces 
%       for finding narrowest voltage peaks
if isempty(vVec3)
    vVec3 = movingaveragefilter(vVec1, smoothWindowMs, siMs);
end

%% Find indices from the time vector
% Index of current application 
idxStimStart = find(tVec0 >= stimStartMs, 1);

% Index to start LTS peak search
idxMinPeakTime = find(tVec0 >= minPeakTimeMs, 1);

% Indices for searching for LTS
%   Note: first index for LTS search must be after current peak
idxSearchBegin = find(tVec0 >= searchWindowMs(1), 1);
idxSearchBegin = max(idxSearchBegin, idxMinPeakTime);
idxSearchEnd = find(tVec0 <= searchWindowMs(2), 1, 'last'); 

% Indices for calculating maxNoise
if computeMaxNoiseFlag
    idxNoiseBegin = find(tVec0 >= noiseWindowMs(1), 1);
    idxNoiseEnd = find(tVec0 <= noiseWindowMs(2), 1, 'last');
    indNoise = idxNoiseBegin:idxNoiseEnd;
end

% Indices for calculating baseline voltage
if computeActVholdFlag
    indBase = idxStimStart - fliplr(1:round(baseWidthMs / siMs));
end

% Print info
if verbose
    fprintf('Time of current application == %g ms\n', stimStartMs);
    fprintf('Index of current application == %d\n', idxStimStart);
    fprintf('Time of minimum LTS peak == %g ms\n', minPeakTimeMs);
    fprintf('Index of minimum LTS peak time == %d\n', idxMinPeakTime);
end

%% Find LTSs and parse features
% Initialization
isSpontaneous = false;      % whether it's a spontaneous spike
isOverridden = false;       % whether the algorithm is overridden
indSpikes = [];        % all spike indices (could be burst or spontaneous spike)
spikesPerPeak = 0;  % all peaks are initially assumed to have no spikes

% Calculate maxNoise
% Note: This is used by legacy:
%       4 * standard deviation of values of the median-filtered then 
%           moving-average-filtered trace is the maximum noise
%       Assuming a Gaussian distribution of noise, should contain 95.45%
if computeMaxNoiseFlag
    maxNoise = 4 * nanstd(vVec3(indNoise));
end
if verbose
    fprintf('Maximum noise == %g mV\n', maxNoise);
end

% Calculate baseline voltage (holding potential)
if computeActVholdFlag
    actVhold = nanmean(vVec1(indBase));
    if verbose
        fprintf('Actual holding potential == %g mV\n', actVhold);
    end
else
    actVhold = NaN;
end

%% Set up the 2nd derivative of median-filtered voltage trace
% Differentiate vVec3, the median-filtered then smoothed voltage trace
dvVec3 = diff(vVec3) / siMs;

% Smooth out dvVec3
dvVec3Sm = movingaveragefilter(dvVec3, smoothWindowMs, siMs);

% Differentiate dvVec3Sm
d2vVec3 = diff(dvVec3Sm) / siMs;

%% Restrict to search window
% Find the indices for the search window
indSearch = idxSearchBegin:idxSearchEnd;
if verbose
    fprintf('Finding peak in the following window: [%g, %g]\n', ...
            siMs*idxSearchBegin, siMs*idxSearchEnd);
end

% Extract voltage vector of interest for detecting LTS & calculating LTS 2ndder
v3 = vVec3(indSearch);               

% Extract voltage vector of interest for detecting spikes (ap)
v0 = vVec0(indSearch);               

% Extract voltage vector of interest for calculating LTS amplitude, 
%   prominence, width
v1 = vVec1(indSearch);               

% Extract 2nd derivative vector of interest
%   Note: two differentiations: 
%           shift right (don't -1) than shift left (-1)
ddv3 = d2vVec3(indSearch - 1);

%% Process all voltage peaks
% TODO: Make a function parse_peaks.m
% Find all local maxima of v3
[peak3Amp, peak3Idx, peak3Width, peak3Prom] = findpeaks(v3);

% Count the local maxima
nPeak3s = numel(peak3Amp);

% If no local maxima exist, find the absolute maximum
if nPeak3s == 0 
    [peak3Amp, peak3Idx] = max(v3);
    peak3Width = 0;
    peak3Prom = 0;
    nPeak3s = 1;
end

% Find peak boundaries and spikes
idxPeak3Start = zeros(nPeak3s, 1);
idxPeak3End = zeros(nPeak3s, 1);
peak3SecDer = zeros(nPeak3s, 1);
peak3Orig = cell(nPeak3s, 1);
indSpikesV0Peak3 = cell(nPeak3s, 1);
idxFirstSpikePeak3 = zeros(nPeak3s, 1);
for iPeak = 1:nPeak3s
    % Get the current peak index
    idxPeak3This = peak3Idx(iPeak);

    % Find peak starting indices
    %   Note: findpeaks() will not work for < 3 points
    if idxPeak3This < 3
        idxPeak3StartThis = 1;
    else
        % Find the first minimum to the left: 
        %   1. Flip and invert 
        %   2. Find first voltage peak
        [~, ind4] = findpeaks(-flipud(v3(1:idxPeak3This)), ...
                            'MinPeakProminence', 0, 'NPeaks', 1);

        % Compute the index relative to v3
        if isempty(ind4)
            idxPeak3StartThis = 1;
        else
            idxPeak3StartThis = (idxPeak3This + 1) - ind4;
        end
    end

    % Find peak ending indices
    %   Note: findpeaks() will not work for < 3 points
    if idxPeak3This > numel(v3) - 2
        idxPeak3EndThis = numel(v3);
    else
        % Find the first minimum to the right: 
        %   1. Flip and invert 
        %   2. Find first voltage peak
        [~, ind5] = findpeaks(-v3(idxPeak3This:end), ...
                            'MinPeakProminence', 0, 'NPeaks', 1);
        if isempty(ind5)
            idxPeak3EndThis = numel(v3);
        else
            idxPeak3EndThis = (idxPeak3This - 1) + ind5;
        end
    end

    % Find the most negative 2nd derivative over the entire peak
    if idxPeak3StartThis == idxPeak3EndThis
        peak3SecDer(iPeak) = ddv3(idxPeak3StartThis);
    else
        peak3SecDer(iPeak) = min(ddv3(idxPeak3StartThis:idxPeak3EndThis));
    end

    % Extract the peak from the original data vector
    peak3Orig{iPeak} = v0(idxPeak3StartThis:idxPeak3EndThis);        

    % Find all "spikes" within the peak
    [spikesAmp, spikesIdx] = findpeaks(peak3Orig{iPeak});

    % Action potentials must be greater than threshold
    stemp1 = find(spikesAmp > spThr);
    if ~isempty(stemp1)
        % Record spike indices relative to v0 or v3
        indSpikesV0Peak3{iPeak} = (idxPeak3StartThis - 1) + spikesIdx(stemp1);

        % Record index of first spike relative to v0 or v3
        idxFirstSpikePeak3(iPeak) = indSpikesV0Peak3{iPeak}(1);
    end

    % Store in arrays
    idxPeak3Start(iPeak) = idxPeak3StartThis;
    idxPeak3End(iPeak) = idxPeak3EndThis;
end

%% Classify each peak and select the LTS peak
% Algorithm for detecting whether there is an LTS and classifying other peaks:
% (1) Eliminate all voltage peaks with prominence < maxNoise
%     then find the first voltage peak with second derivative below threshold.
%     If it has spikes, 
%       it is an LTS only if the first spike occurs before the peak.
% (2) If doesn't exist, 
%       eliminate all voltage peaks with prominence < maxNoise, 
%       then find the narrowest voltage peak.
% (3) If still doesn't exist, find the narrowest voltage peak.

% Find all peaks with prominence greater than maximum noise
%   Note: previously 1 mV, now different for each trace
ptemp1 = find(peak3Prom > maxNoise);

% Find all peaks with the second derivative reaching (<=) threshold
ptemp2 = find(peak3SecDer <= ltsThr);

% Find all peaks that are LTS candidates by threshold
ptemp3 = intersect(ptemp1, ptemp2);

if isempty(ptemp1)        % Condition (3)
    % Not LTS by prominence
    ltsPeakValue = NaN;
    ltsPeakTime = NaN;
    maxSlopeTime = NaN;
    maxSlopeValue = NaN;

    % Find the narrowest voltage peak
    [peak2ndDer, iPeak3Sel] = min(peak3SecDer);
    peakClass = 1;
    peakClassLabel = ['not LTS: peak prominence ', ...
                    num2str(peak3Prom(iPeak3Sel)), ...
                    ' mV <= maximum noise ', ...
                    num2str(maxNoise), ' mV'];
    
    % Override the algorithm if 3/4 experts think it's an LTS
    if ismember(fileBase, fileBasesToOverride{idxLooksLikeLtsNotByProminence})
        isOverridden = true;
        peakClass = 5;
        peakClassLabel = ['LTS by overrule: peak prominence ', ...
                        num2str(peak3Prom(iPeak3Sel)), ...
                        ' mV <= maximum noise ', ...
                        num2str(maxNoise), ' mV'];
    end
elseif isempty(ptemp3) ...    % Condition (2)
    % Not LTS by narrowness
    ltsPeakValue = NaN;
    ltsPeakTime = NaN;
    maxSlopeTime = NaN;
    maxSlopeValue = NaN;

    % Find the narrowest voltage peak with prominence greater than maximum noise
    [peak2ndDer, pp1] = min(peak3SecDer(ptemp1));
    iPeak3Sel = ptemp1(pp1);
    peakClass = 2;
    peakClassLabel = ['not LTS: 2nd derivative ', ...
                        num2str(peak2ndDer), ...
                        ' V^2/s^2 > LTS threshold ', ...
                        num2str(ltsThr), ' V^2/s^2'];

    % Override the algorithm if 3/4 of experts think it's an LTS
    if ismember(fileBase, fileBasesToOverride{idxLooksLikeLtsNotByNarrowness})
        isOverridden = true;
        peakClass = 5;
        peakClassLabel = ['LTS by overrule: 2nd derivative ', ...
                            num2str(peak2ndDer), ...
                            ' V^2/s^2 > LTS threshold ', ...
                            num2str(ltsThr), ' V^2/s^2'];
    end

    % There were a few cases (see 'Looks_like_missed_LTS') 
    %   where the correct LTS was missed by prominence
    %   The first narrow peak regardless of prominence 
    %       is the correct one in this case
    if ismember(fileBase, fileBasesToOverride{idxLooksLikeMissedLts})
        iPeak3Sel = ptemp2(1);
        peak2ndDer = peak3SecDer(iPeak3Sel);
        isOverridden = true;
        peakClass = 5;
        peakClassLabel = ['LTS by overrule: peak prominence ', ...
                            num2str(peak3Prom(iPeak3Sel)), ...
                            ' mV <= maximum noise ', ...
                            num2str(maxNoise), ' mV'];
    end
else                        % Condition (1)
    % Select the first peak that is an LTS candidate by prominence & 2nd der
    %   Note: There is one case (see 'Missed_LTS_by_order', F092710_0006_25) 
    %           where the second peak is the correct LTS
    if ismember(fileBase, fileBasesToOverride{idxMissedLtsByOrder})
        iPeak3Sel = ptemp3(2);
    else
        iPeak3Sel = ptemp3(1);
    end
    peak2ndDer = peak3SecDer(iPeak3Sel);

    % Check whether it's a spontaneous spike
    %     (Based on following observation of shape:
    %         LTS:         first spike occurs before "LTS" peak on mfmaf trace, 
    %            except in four cases (see 'Missed_LTS_by_shape'),
    %             where the first spike occurred after the peak
    %     spontaneous spike:     first spike occurs after "LTS" peak on mfmaf trace)
    %   Note: minSp2PkTime is obsolete: currently set to 0
    sp2pkSamples = round(minSp2PkTime/siMs);
    if ~ismember(fileBase, fileBasesToOverride{idxMissedLtsByShape}) ...
            && idxFirstSpikePeak3(iPeak3Sel) ~= 0 ...
            && peak3Idx(iPeak3Sel) < idxFirstSpikePeak3(iPeak3Sel) + sp2pkSamples
        % Not LTS by shape
        ltsPeakValue = NaN;
        ltsPeakTime = NaN;
        maxSlopeTime = NaN;
        maxSlopeValue = NaN;

        isSpontaneous = false;
        peakClass = 3;
        peakClassLabel = ['not LTS: peak index ', ...
                            num2str(peak3Idx(iPeak3Sel)), ...
                            ' < index of first spike ', ...
                            num2str(idxFirstSpikePeak3(iPeak3Sel))];
    end
end
if verbose
    fprintf('Selected peak 2nd derivative == %g V^2/s^2\n', peak2ndDer);
end

%% Analyze the selected LTS peak
% Find all local maxima
[peak1Amp, peak1Idx, peak1Width, peak1Prom] = findpeaks(v1); 

% Count local maxima
nPeak1s = numel(peak1Amp);

% If no local maxima exist, find the absolute maximum
if nPeak1s == 0
    [peak1Amp, peak1Idx] = max(v1);
    peak1Width = 0;
    peak1Prom = 0;
    nPeak1s = 1;
end

% Get the index of the selected peak in v0, v1 or v3
idxPeakV3 = peak3Idx(iPeak3Sel);

% Find the peak of v1 that corresponds to the selected peak of v3
[~, iPeak1Sel] = min(abs(peak1Idx - idxPeakV3));

% Record the prominence of selected peak
% this is the "minimum vertical distance that the signal must 
% descend on either side of the peak before either climbing back 
% to a level higher than the peak or reaching an endpoint" in mV
peakProm = peak1Prom(iPeak1Sel);

% Record the width of the peak at half-prominence in ms
peakWidth = peak1Width(iPeak1Sel) * siMs;

% Record LTS peak index relative to vVec0 or vVec1
idxPeak = (idxSearchBegin - 1) + peak1Idx(iPeak1Sel);

% Compute LTS peak time (delay) in ms
peakTime = (idxPeak - idxStimStart) * siMs;

% Print messages
if verbose
    fprintf('Selected peak prominence == %g mV\n', peakProm);
    fprintf('Selected peak width == %g ms\n', peakWidth);
    fprintf('Selected peak time == %g ms\n', peakTime);
end

% Compute LTS peak end points
%   Note: The following may be changed later for bursts
idxPeakStart = (idxSearchBegin - 1) + idxPeak3Start(iPeak3Sel);
idxPeakEnd = (idxSearchBegin - 1) + idxPeak3End(iPeak3Sel);

% Save old peak end points for plotting
idxPeakStartOld = idxPeakStart;
idxPeakEndOld = idxPeakEnd;

% Initialize spike indices with what was found before
indSpikesV0 = indSpikesV0Peak3{iPeak3Sel};

% Record spike indices relative to vVec0 & count spikes per peak
if ~isempty(indSpikesV0)
    % Get all spike indices relative to vVec0
    indSpikes = (idxSearchBegin - 1) + indSpikesV0;

    % Count spikes per peak
    spikesPerPeak = numel(indSpikesV0);
end

% Find LTS amplitude, delay, maximum slope delay and value
%   if either peakClass hasn't been classified yet or it's an 'LTS by overrule'
if (~isempty(ptemp3) && ~isSpontaneous) || isOverridden
    % Compute LTS amplitude, using the median-filtered voltage trace
    ltsPeakValue = peak1Amp(iPeak1Sel);

    % Compute LTS delay, using median-filtered voltage trace
    ltsPeakTime = peakTime;

    % Compute threshold for spike detection (mV)
    spikeThreshold = ltsPeakValue + spThrRelLts;

    % Find approximate max slope based off of v3
    spacingSize = round(slopeSpacing/siMs);                           % number of indices apart to calculate slope
    v3pkLeft = v3(idxPeak3Start(iPeak3Sel):idxPeak3End(iPeak3Sel) - spacingSize);     % the left points for slope calculation
    v3pkRight = v3(idxPeak3Start(iPeak3Sel) + spacingSize:idxPeak3End(iPeak3Sel));    % the right points for slope calculation
    v3pkSlopes = (v3pkRight - v3pkLeft)/slopeSpacing;               % all slopes in V/s
    [maxSlopeValAppr, ~] = max(v3pkSlopes);                           % find approximate maximum slope

    % Moving-average-filter median-filtered traces for calculating maximum slope
    mafw3 = mafw3Dv/maxSlopeValAppr;      % width in ms for the moving average filter for finding slopes
    v4 = movingaveragefilter(v1, mafw3, siMs);  % voltage vector of interest for calculating maxslope

    % Find max slope based off of v4
    spacingSize = round(slopeSpacing/siMs);                       % number of indices apart to calculate slope
    v4pkLeft = v4(idxPeak3Start(iPeak3Sel):idxPeak3End(iPeak3Sel) - spacingSize); % the left points for slope calculation
    v4pkRight = v4(idxPeak3Start(iPeak3Sel) + spacingSize:idxPeak3End(iPeak3Sel));% the right points for slope calculation
    v4pkSlopes = (v4pkRight - v4pkLeft)/slopeSpacing;           % all slopes in V/s
    [maxSlopeValue, tempInd] = max(v4pkSlopes);                     % find maximum slope

    maxslopeind = tempInd + spacingSize + (idxPeakStart - 1);           % the maxslope index is the right point
    maxSlopeTime = (maxslopeind - idxStimStart) * siMs;               % delay in ms of maximum slope after IPSC starts
    peakFeatureLabel = ['max slope = ', num2str(maxSlopeValue), ' V/s; ', ...
                        'prominence = ', num2str(peakProm), ' mV; ', ...
                        'width = ', num2str(peakWidth), ' ms'];
end

% Determine whether it's a "burst" and count action potentials (spikes)
if isnan(ltsPeakTime)           % must be an "LTS" to begin with
    burstTime = NaN;
    spikesPerBurst = NaN;
    spikeThreshold = NaN;
    firstSpikeTime = NaN;
    lastSpikeTime = NaN;
    maxSpikeAmp = NaN;
    minSpikeAmp = NaN;
    spikeFrequency = NaN;
    spikeAdaptation = NaN;
else
    if isempty(indSpikes)          % no spikes detected, not a "burst"
        burstTime = NaN;
        spikesPerBurst = NaN;
        spikeThreshold = NaN;
        firstSpikeTime = NaN;
        lastSpikeTime = NaN;
        maxSpikeAmp = NaN;
        minSpikeAmp = NaN;
        spikeFrequency = NaN;
        spikeAdaptation = NaN;
        if peak2ndDer <= ltsThrAlt && ~isOverridden
            peakClass = 8;
            peakClassLabel = ['LTS with no burst; ', ...
                        'definite: 2nd der ', ...
                        num2str(peak2ndDer), ...
                        ' V^2/s^2'];
        elseif ~isOverridden
            peakClass = 6;
            peakClassLabel = ['LTS with no burst; ', ...
                        'contentious: 2nd der ', ...
                        num2str(peak2ndDer), ...
                        ' V^2/s^2 > ', ...
                        'alt LTS thr ', ...
                        num2str(ltsThrAlt), ...
                        ' V^2/s^2'];
        end
    else
        if peak2ndDer <= ltsThrAlt && ~isOverridden
            peakClass = 9;
            peakClassLabel = ['LTS with burst; definite: ', ...
                        '2nd der ', ...
                        num2str(peak2ndDer), ...
                        ' V^2/s^2'];
        elseif ~isOverridden
            peakClass = 7;
            peakClassLabel = ['LTS with burst; contentious: ', ...
                        '2nd der ', ...
                        num2str(peak2ndDer), ...
                        ' V^2/s^2 > ', ...
                        'alt LTS thr ', ...
                        num2str(ltsThrAlt), ...
                        ' V^2/s^2'];
        end

        if ~ismember(fileBase, fileBasesToOverride{idxSpikesPerBurstIncorrect}) ...
            % Re-detect spikes by redefining peak bounds using maxNoise as MinPeakProminence
            % Update peak starting indices
            if idxPeakV3 < 3                  % findpeaks will not work for < 3 points
                idxPeak3Start(iPeak3Sel) = 1;
            else
                [~, ind4] = findpeaks(-flipud(v3(1:idxPeakV3)), ...
                        'MinPeakProminence', maxNoise, 'NPeaks', 1);
                    % first minimum to the left: flip and invert, 
                    % then find first voltage peak with prominence 
                    % >= maxNoise (previously 2 mV, now different for each trace)
                if isempty(ind4)
                    idxPeak3Start(iPeak3Sel) = 1;
                else
                    idxPeak3Start(iPeak3Sel) = (idxPeakV3 + 1) - ind4;
                end
            end
            % Update peak ending indices
            if idxPeakV3 > numel(v3) - 2    % findpeaks will not work for < 3 points
                idxPeak3End(iPeak3Sel) = numel(v3);
            else
                [~, ind5] = findpeaks(-v3(idxPeakV3:end), ...
                        'MinPeakProminence', maxNoise, 'NPeaks', 1);
                    % first minimum to the right: invert, 
                    % then find first voltage peak with prominence 
                    % >= maxNoise (previously 2 mV, now different for each trace)
                if isempty(ind5)
                    idxPeak3End(iPeak3Sel) = numel(v3);
                else
                    idxPeak3End(iPeak3Sel) = (idxPeakV3 - 1) + ind5;
                end
            end
            % Detect spikes in original trace within the peak
            peak3Orig{iPeak3Sel} = v0(idxPeak3Start(iPeak3Sel):idxPeak3End(iPeak3Sel)); % take only the peak part
            [spikesAmp, spikesIdx] = findpeaks(peak3Orig{iPeak3Sel});   % find all "spikes" within the peak
            stemp1 = find(spikesAmp > spikeThreshold);          % action potentials must be greater than actual threshold
                                % Note: this is a different threshold than before,
                                %    so could potentially change the classification of the peak

            % Update spike indices relative to v0 or v3
            if ~isempty(stemp1)
                indSpikesV0 = (idxPeak3Start(iPeak3Sel) - 1) + spikesIdx(stemp1);
            else
                indSpikesV0 = [];
            end

            % Update peak bound indices relative to vVec0
            idxPeakStart = (idxSearchBegin - 1) + idxPeak3Start(iPeak3Sel);       % narrowest peak lower bound index
            idxPeakEnd = (idxSearchBegin - 1) + idxPeak3End(iPeak3Sel);       % narrowest peak upper bound index
        end

        % Two cases: spikes are still found or not
        if ~isempty(indSpikesV0)
            % Update spike indices and spikes per peak
            indSpikes = (idxSearchBegin - 1) + indSpikesV0;            % spike indices
            spikesPerPeak = numel(indSpikes);                 % count spikes per peak

            % Record time of first and last spike
            firstSpikeTime = (indSpikes(1) - idxStimStart) * siMs;
            lastSpikeTime = (indSpikes(end) - idxStimStart) * siMs;
            spikeTimeLabel = ['spikes; first = ', num2str(firstSpikeTime), ...
                                ' ms; last = ', num2str(lastSpikeTime), ' ms'];

            % Spikes per burst is the the same spikes per peak but not zero
            spikesPerBurst = spikesPerPeak;                 % spikes per burst

            % Define the burst region
            vBurst0 = peak3Orig{iPeak3Sel};                         % voltage trace of burst
            idxFirstSpikePeak = idxFirstSpikePeak3(iPeak3Sel) - idxPeak3Start(iPeak3Sel) + 1;% index of first spike in vBurst0

            % Find first minimum to the left of first spike
            if idxFirstSpikePeak < 3                       % findpeaks will not work for < 3 points
                idxFirstSpikeStart = 1;
                amp6 = -vBurst0(1);
            else
                [amp6, ind6] = findpeaks(-flipud(vBurst0(1:idxFirstSpikePeak)), 'NPeaks', 1);
                if isempty(ind6)
                    idxFirstSpikeStart = 1;
                    amp6 = -vBurst0(1);
                else
                    idxFirstSpikeStart = (idxFirstSpikePeak + 1) - ind6;
                end
            end
            % Find first minimum to the right of first spike
            if idxFirstSpikePeak > numel(vBurst0) - 2     % findpeaks will not work for < 3 points
                idxFirstSpikeEnd = numel(vBurst0);
                amp7 = -vBurst0(end);
            else
                [amp7, ind7] = findpeaks(-vBurst0(idxFirstSpikePeak:end), 'NPeaks', 1);
                if isempty(ind7)
                    idxFirstSpikeEnd = numel(vBurst0);
                    amp7 = -vBurst0(end);
                else
                    idxFirstSpikeEnd = (idxFirstSpikePeak - 1) + ind7;
                end
            end
            amp8 = max([-amp6 -amp7]);          % take the higher of the two minimums

            % Find the last index lower than the base of the first spike in terms of vBurst0
            ind8 = find(vBurst0(1:idxFirstSpikePeak) < amp8, 1, 'last');

            % Find the burst onset time (delay)
            if isempty(ind8) || ind8 < 4
                % Burst onset index is the beginning of the LTS peak in terms of vVec0
                idxBurstOnset = (idxSearchBegin - 1) + idxPeak3Start(iPeak3Sel);
            else                        % trace to differentiate must have at least 4 elements
                % Find the corresponding median-filtered then smoothed trace for the burst region
                vBurst4 = v4(idxPeak3Start(iPeak3Sel):idxPeak3End(iPeak3Sel));
                
                % Find the index (in terms of vBurst4) of the maximum of the 3rd derivative of
                %   the voltage trace in between the start of burst and the last point before the first spike
                %   Add 3 indices to account for the loss of an element on the left after each diff()
                [~, ind9] = max(diff(diff(diff(vBurst4(1:ind8)))) / (siMs^3));
                max3rdderi = ind9 + 3;
                
                % Burst onset index is in terms of vVec0
                idxBurstOnset = (idxSearchBegin - 1) + (idxPeak3Start(iPeak3Sel) - 1) + max3rdderi;
            end
            burstTime = (idxBurstOnset - idxStimStart) * siMs;    % burst onset time (delay)
            burstTimeLabel = ['burst onset; delay = ', num2str(burstTime), ' ms'];

            % The computed spike threshold is the voltage at burst onset
            spikeThreshold = vVec0(idxBurstOnset);

            % Find the maximum spike amplitude, minimum spike amplitude, 
            %   spike frequency and spike adaptation
            if numel(indSpikes) >= 1
                maxSpikeAmp = max(vVec0(indSpikes));
                minSpikeAmp = min(vVec0(indSpikes));
                if numel(indSpikes) >= 2
                    % Compute the spike frequency in Hz
                    spikeFrequency = compute_spike_frequency(indSpikes, siMs);

                    if numel(indSpikes) >= 3
                        spikeAdaptation = 100 * (tVec0(indSpikes(end)) - tVec0(indSpikes(end-1))) ...
                                            / (tVec0(indSpikes(2)) - tVec0(indSpikes(1)));
                                    % spike adaptation (%) is last ISI/first ISI
                    else
                        spikeAdaptation = NaN;
                    end
                else
                    spikeFrequency = NaN;
                    spikeAdaptation = NaN;
                end
            else
                maxSpikeAmp = NaN;
                minSpikeAmp = NaN;
                spikeFrequency = NaN;
                spikeAdaptation = NaN;
            end
        else                        % no longer a burst
            indSpikes = [];
            spikesPerPeak = 0;
            burstTime = NaN;
            spikesPerBurst = NaN;
            spikeThreshold = NaN;
            firstSpikeTime = NaN;
            lastSpikeTime = NaN;
            maxSpikeAmp = NaN;
            minSpikeAmp = NaN;
            spikeFrequency = NaN;
            spikeAdaptation = NaN;
            if peak2ndDer <= ltsThrAlt && ~isOverridden
                peakClass = 8;
                peakClassLabel = ['LTS with no burst; ', ...
                            'definite: 2nd der ', ...
                            num2str(peak2ndDer), ...
                            ' V^2/s^2'];
            elseif ~isOverridden
                peakClass = 6;
                peakClassLabel = ['LTS with no burst; ', ...
                            'contentious: 2nd der ', ...
                            num2str(peak2ndDer), ...
                            ' V^2/s^2 > ', ...
                            'alt LTS thr ', ...
                            num2str(ltsThrAlt), ...
                            ' V^2/s^2'];
            end
        end
    end
end

%% Deal with false positives
if ismember(fileBase, fileBasesToOverride{idxNoiseInTrace}) ...
    || ismember(fileBase, fileBasesToOverride{idxSpontLtsOrBurst}) ...
    || ismember(fileBase, fileBasesToOverride{idxWideLtsCouldBeNoise})
    % Not LTS by overrule
    ltsPeakValue = NaN;
    ltsPeakTime = NaN;
    maxSlopeTime = NaN;
    maxSlopeValue = NaN;
    burstTime = NaN;
    spikesPerBurst = NaN;
    spikeThreshold = NaN;
    firstSpikeTime = NaN;
    lastSpikeTime = NaN;
    maxSpikeAmp = NaN;
    minSpikeAmp = NaN;
    spikeFrequency = NaN;
    spikeAdaptation = NaN;
    isSpontaneous = false;

    % Label info according to the previously detected peakClass
    if peakClass == 6 || peakClass == 7
        peakClassLabel = ['not LTS by overrule: 2nd der ', ...
                            num2str(peak2ndDer), ' V^2/s^2 < ', ...
                            'LTS threshold ', num2str(ltsThr), ' V^2/s^2'];
    elseif peakClass == 8 || peakClass == 9
        peakClassLabel = ['not LTS by overrule: 2nd der ', ...
                            num2str(peak2ndDer), ' V^2/s^2 < ', ...
                            'alt LTS thr ', num2str(ltsThrAlt), ' V^2/s^2'];
    end
    
    % Make this peakClass == 4
    peakClass = 4;            
end

% Decide whether LTS could have been missed
if isempty(ptemp3) && ~isempty(ptemp2)
    couldHaveMissed = true;
else
    couldHaveMissed = false;
end

% fprintf('Selected peak class == %d\n', peakClass);
% fprintf('Spikes per peak == %d\n', spikesPerPeak);
% fprintf('LTS peak time == %g ms\n', ltsPeakTime);
% fprintf('LTS amplitude (absolute) == %g mV\n', ltsPeakValue);
% fprintf('LTS maximum slope time == %g ms\n', maxSlopeTime);
% fprintf('LTS maximum slope amplitude == %g V/s\n', maxSlopeValue);
% fprintf('Burst onset time == %g ms\n', burstTime);
% fprintf('Spikes per burst == %d\n\n', spikesPerBurst);

% Store in parsedParams structure
parsedParams.fileBase = fileBase;
parsedParams.siMs = siMs;
parsedParams.stimStartMs = stimStartMs;
parsedParams.minPeakTimeMs = minPeakTimeMs;
parsedParams.medfiltWindowMs = medfiltWindowMs;
parsedParams.smoothWindowMs = smoothWindowMs;
parsedParams.baseWidthMs = baseWidthMs;
parsedParams.ltsThr = ltsThr;
parsedParams.ltsThrAlt = ltsThrAlt;
parsedParams.spThr = spThr;
parsedParams.spThrRelLts = spThrRelLts;
parsedParams.siMsRes = siMsRes;
parsedParams.minSp2PkTime = minSp2PkTime;
parsedParams.slopeSpacing = slopeSpacing;
parsedParams.mafw3Dv = mafw3Dv;
parsedParams.slopeSegYHalf = slopeSegYHalf;
parsedParams.actVhold = actVhold;
parsedParams.maxNoise = maxNoise;
parsedParams.peakTime = peakTime;
parsedParams.peak2ndDer = peak2ndDer;
parsedParams.peakProm = peakProm;
parsedParams.peakWidth = peakWidth;
parsedParams.peakClass = peakClass;
parsedParams.spikesPerPeak = spikesPerPeak;
parsedParams.ltsPeakTime = ltsPeakTime;
parsedParams.ltsPeakValue = ltsPeakValue;
parsedParams.maxSlopeTime = maxSlopeTime;
parsedParams.maxSlopeValue = maxSlopeValue;
parsedParams.burstTime = burstTime;
parsedParams.spikesPerBurst = spikesPerBurst;
parsedParams.spikeThreshold = spikeThreshold;
parsedParams.firstSpikeTime = firstSpikeTime;
parsedParams.lastSpikeTime = lastSpikeTime;
parsedParams.maxSpikeAmp = maxSpikeAmp;
parsedParams.minSpikeAmp = minSpikeAmp;
parsedParams.spikeFrequency = spikeFrequency;
parsedParams.spikeAdaptation = spikeAdaptation;
parsedParams.isSpontaneous = isSpontaneous;
parsedParams.isOverridden = isOverridden;
parsedParams.couldHaveMissed = couldHaveMissed;

% Store in parsedData structure
parsedData.tVec0 = tVec0;
parsedData.tVec2 = tVec2;
parsedData.vVec0 = vVec0;
parsedData.vVec1 = vVec1;
parsedData.vVec2 = vVec2;
parsedData.vVec3 = vVec3;
parsedData.noiseWindowMs = noiseWindowMs;
parsedData.searchWindowMs = searchWindowMs;
if computeActVholdFlag
    parsedData.indBase = indBase;
end
parsedData.indSearch = indSearch;

%% Plot figures
if plotFlag
    % Turn off legend autoupdate for back compatibility
    defaultLegendAutoUpdateOrig = get(groot, 'defaultLegendAutoUpdate');
    set(groot, 'defaultLegendAutoUpdate', 'off');

    % Check if needed outSubDirs exist in outFolder
    check_subdir(outFolder, outSubDirs);

    % Get the end points of the baseline detection
    if computeActVholdFlag
        idxBaseStart = indBase(1);
        idxBaseEnd = indBase(end);
    end

    % Create a file base string for titles
    fileBaseTitle = replace(fileBase, '_', '\_');

    % Compute info needed for plotting peak properties
    if ~isnan(ltsPeakTime)
        % Compute maxslope line segment limits
%{
        leftslopepoint_y = min(vVec3(idxPeakEnd), vVec3(maxslopeind) - slopeSegYHalf);
        rightslopepoint_y = min(vVec3(idxPeak), vVec3(maxslopeind) + slopeSegYHalf);
        leftslopepoint_x = (leftslopepoint_y - vVec3(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);
        rightslopepoint_x = (rightslopepoint_y - vVec3(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);
%}

        leftslopepoint_y = min(vVec1(idxPeakEnd), vVec1(maxslopeind) - slopeSegYHalf);
        rightslopepoint_y = min(vVec1(idxPeak), vVec1(maxslopeind) + slopeSegYHalf);
        leftslopepoint_x = (leftslopepoint_y - vVec1(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);
        rightslopepoint_x = (rightslopepoint_y - vVec1(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);

        % Compute info for peak prominence line segment
        pkprom_x = tVec0(idxPeak);            % time value (ms) of peak prominence
        pkprom_y1 = vVec1(idxPeak) - peakProm;    % voltage value (mV) of bottom of peak
        pkprom_y2 = vVec1(idxPeak);            % voltage value (mV) of top of peak

        % Compute info for peak width line segment
        halfprom_y = vVec1(idxPeak)- peakProm/2;            % voltage value (mV) at half prominence
        [~, tempInd] = min(abs(vVec1(idxPeakStart:idxPeak) - halfprom_y));    % Find the index on the rising phase 
                                        % of the peak whose voltage value 
                                        % is closest to half prominence
        pkwidth_leftind = tempInd + idxPeakStart - 1;        % peak width segment left index relative to vVec1
        pkwidth_x1 = tVec0(pkwidth_leftind);            % time value (ms) of left end
        pkwidth_x2 = tVec0(pkwidth_leftind) + peakWidth;    % time value (ms) of right end
    end

    % Plot voltage traces
    % fprintf('Plotting voltage trace ...\n');
    h = figure(5000);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Voltage trace');
    clf(h);
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
    plot(tVec0, vVec0, 'b-', 'LineWidth', 0.5); hold on;
    plot(tVec0, vVec1, 'g-', 'LineWidth', 0.5);
    plot(tVec0, vVec3, 'r-', 'LineWidth', 0.5);
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)    % LTS with bursts
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), vVec1(idxPeak), 'go', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'ro', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                % line for maxslope
        plot(tVec0(idxBurstOnset), vVec0(idxBurstOnset), 'g>', 'MarkerSize', 10);       % triangle for burst onset time
        plot(tVec0(indSpikes), vVec0(indSpikes), 'gx', 'MarkerSize', 10);     % crosses for spikes
        legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakFeatureLabel, burstTimeLabel, spikeTimeLabel, ...
            'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
        line([tVec0(1), tVec0(idxBurstOnset)], [spikeThreshold, spikeThreshold], ...
                'LineStyle', ':');                                      % line for spike threshold
        line([tVec0(1), tVec0(end)], [maxSpikeAmp, maxSpikeAmp], 'LineStyle', '--');    % line for maxSpikeAmp
        line([tVec0(1), tVec0(end)], [minSpikeAmp, minSpikeAmp], 'LineStyle', '--');    % line for minSpikeAmp
        text(.05, .15, ['Spike Frequency: ', num2str(spikeFrequency), ' Hz'], 'Units', 'normalized');
        text(.05, .1, ['Spike Adaptation: ',  num2str(spikeAdaptation), '%'], 'Units', 'normalized');
        text(.05, .05, ['Spike Threshold: ',  num2str(spikeThreshold), 'mV'], 'Units', 'normalized');
    elseif ~isnan(ltsPeakTime)            % LTS without bursts
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), vVec1(idxPeak), 'bo', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(idxPeak), vVec1(idxPeak), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                             % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                              % line for maxslope
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
        legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakFeatureLabel, 'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
    else                        % not LTS
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(idxPeak), vVec1(idxPeak), 'kx', 'MarkerSize', 10);
        end
        if isempty(indSpikes)            % noise
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'Location', 'SouthOutside')
        elseif isSpontaneous            % spontaneous spikes
            plot(tVec0(indSpikes(1)), vVec0(indSpikes(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike', ...
                'Location', 'SouthOutside')
        else                    % spontaneous spike or noise
            plot(tVec0(indSpikes(1)), vVec0(indSpikes(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike or noise', ...
                'Location', 'SouthOutside')
        end
    end
    if computeActVholdFlag
        plot(tVec0(idxBaseStart), vVec1(idxBaseStart), 'g>');
        plot(tVec0(idxBaseEnd), vVec1(idxBaseEnd), 'y<');
    end
    plot(tVec0(idxPeakStartOld), vVec3(idxPeakStartOld), 'k*');
    plot(tVec0(idxPeakEndOld), vVec3(idxPeakEndOld), 'k*');
    plot(tVec0(idxPeakStart), vVec3(idxPeakStart), 'r*');
    plot(tVec0(idxPeakEnd), vVec3(idxPeakEnd), 'r*');
    xlim([tVec0(1) tVec0(end)]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(fileBaseTitle);
    figname = fullfile(outFolder, outSubDirs{1}, [fileBase, '.png']);
    saveas(h, figname);

    % Plot scaled voltage trace for better comparison
    figure(h);
%    ylim([-120 20]);
    plot(tVec0, vVec1, 'g-', 'LineWidth', 1);
    xlimits = searchWindowMs;
    xlim(xlimits);
    ylim([-100 -40]);            % Fix y-axis to determine whether the trace is good for fitting
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
    legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', peakClassLabel);
    figname = fullfile(outFolder, outSubDirs{4}, [fileBase, '_scaled.png']);
    saveas(h, figname);
    if peak2ndDer > ltsThrAlt ...
        && peak2ndDer <= ltsThr    % in "gray area"
        figname = fullfile(outFolder, outSubDirs{5}, [fileBase, '_scaled.png']);
        saveas(h, figname);
    end
    if couldHaveMissed
        figname = fullfile(outFolder, outSubDirs{6}, [fileBase, '_scaled.png']);
        saveas(h, figname);
    end


    hold off;
%    close(h);

    % Plot LTS analysis
    % fprintf('Plotting LTS analysis ...\n');
    xlimits = searchWindowMs;
    h = figure(5001);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'LTS analysis, moving-average-filtered trace');
    clf(h);
    subplot(3,1,1) % voltage trace
    plot(tVec0, vVec3, 'r-', 'LineWidth', 0.5); hold on
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), vVec1(idxPeak), 'go', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'gd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(idxPeak), vVec1(idxPeak), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'md', 'MarkerSize', 8);
        end
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), vVec1(idxPeak), 'bo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'bd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(idxPeak), vVec1(idxPeak), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'md', 'MarkerSize', 8);
        end
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(idxPeak), vVec1(idxPeak), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tVec0(idxPeakStart), vVec3(idxPeakStart), 'r*');
    plot(tVec0(idxPeakEnd), vVec3(idxPeakEnd), 'r*');
    plot(tVec0(idxPeakStartOld), vVec3(idxPeakStartOld), 'k*');
    plot(tVec0(idxPeakEndOld), vVec3(idxPeakEndOld), 'k*');
    xlim(xlimits);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(['LTS analysis for ', fileBaseTitle, ', moving-average-filtered trace']);
    subplot(3,1,2) % 1st derivative of voltage trace
    plot(tVec0(2:end), dvVec3, 'k-', 'LineWidth', 0.5); hold on
    plot(tVec0(2:end), dvVec3Sm, 'r-', 'LineWidth', 0.5);
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'go', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvVec3(maxslopeind), 'gd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvVec3(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvVec3(maxslopeind), 'md', 'MarkerSize', 8);
        end
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'bo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvVec3(maxslopeind), 'bd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvVec3(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvVec3(maxslopeind), 'md', 'MarkerSize', 8);
        end
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(idxPeak), dvVec3(idxPeak), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tVec0(idxPeakStart), dvVec3(idxPeakStart), 'r*');
    plot(tVec0(idxPeakEnd), dvVec3(idxPeakEnd), 'r*');
    plot(tVec0(idxPeakStartOld), dvVec3(idxPeakStartOld), 'k*');
    plot(tVec0(idxPeakEndOld), dvVec3(idxPeakEndOld), 'k*');
    xlim(xlimits);
    xlabel('Time (ms)');
    ylabel('dV/dT');
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
    legend('unsmoothed', 'smoothed');
    subplot(3,1,3) % 2nd derivative of voltage trace
    plot(tVec0(3:end), d2vVec3, 'k-', 'LineWidth', 0.5); hold on
    line(xlimits, [ltsThr ltsThr], 'Color', 'r', ...
        'LineStyle', '--', 'LineWidth', 0.5);                       % mark LTS threshold
    line(xlimits, [ltsThrAlt ltsThrAlt], 'Color', 'b', ...
        'LineStyle', '--', 'LineWidth', 0.5);                       % mark alt LTS threshold
    xlim(xlimits);
    ax = gca;
    ylimits = get(ax, 'YLim');
    line([tVec0(idxSearchBegin) tVec0(idxSearchBegin)], ...
        ylimits, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);% mark start of search window
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'go', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr        % in "gray area"
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr         % LTS by overrule
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'mo', 'MarkerSize', 10);
        end
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'bo', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr        % in "gray area"
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr         % LTS by overrule
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'mo', 'MarkerSize', 10);
        end
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr        % in "gray area"
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(idxPeak), d2vVec3(idxPeak), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tVec0(idxPeakStart), d2vVec3(idxPeakStart), 'r*');
    plot(tVec0(idxPeakEnd), d2vVec3(idxPeakEnd), 'r*');
    plot(tVec0(idxPeakStartOld), d2vVec3(idxPeakStartOld), 'k*');
    plot(tVec0(idxPeakEndOld), d2vVec3(idxPeakEndOld), 'k*');
    xlabel('Time (ms)');
    ylabel('d2V/dT2');
    figname = fullfile(outFolder, outSubDirs{2}, [fileBase, '_LTSanalysis', '.png']);
    saveas(h, figname);
    hold off;
    % close(h);

    % Plot burst analysis
    % fprintf('Plotting burst analysis ...\n');
    h = figure(5002);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Burst analysis, original trace');
    clf(h);
    plot(tVec0(idxPeakStart:idxPeakEnd), vVec0(idxPeakStart:idxPeakEnd), 'b-', 'LineWidth', 0.5); hold on
    plot(tVec2, vVec2, 'g.', 'LineWidth', 0.5);
    plot(tVec0(idxPeakStart:idxPeakEnd), vVec3(idxPeakStart:idxPeakEnd), 'r-', 'LineWidth', 0.5);
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), vVec1(idxPeak), 'go', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(idxPeak), vVec1(idxPeak), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                             % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                              % line for maxslope
        plot(tVec0(idxBurstOnset), vVec0(idxBurstOnset), 'g>', 'MarkerSize', 10);       % triangle for burst onset time
        plot(tVec0(indSpikes), vVec0(indSpikes), 'gx', 'MarkerSize', 10);     % crosses for spikes
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
        legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakFeatureLabel, burstTimeLabel, spikeTimeLabel, ...
            'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
        line([tVec0(1), tVec0(idxBurstOnset)], [spikeThreshold, spikeThreshold], ...
                'LineStyle', ':');                                      % line for spike threshold
        line([tVec0(1), tVec0(end)], [maxSpikeAmp, maxSpikeAmp], 'LineStyle', '--');    % line for maxSpikeAmp
        line([tVec0(1), tVec0(end)], [minSpikeAmp, minSpikeAmp], 'LineStyle', '--');    % line for minSpikeAmp
        text(.05, .7, ['Spike Frequency: ', num2str(spikeFrequency), ' Hz'], 'Units', 'normalized');
        text(.05, .65, ['Spike Adaptation: ',  num2str(spikeAdaptation), '%'], 'Units', 'normalized');
        text(.05, .6, ['Spike Threshold: ',  num2str(spikeThreshold), 'mV'], 'Units', 'normalized');
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(idxPeak), vVec1(idxPeak), 'bo', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(idxPeak), vVec1(idxPeak), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                % line for maxslope
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
        legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakFeatureLabel, 'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(idxPeak), vVec1(idxPeak), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(idxPeak), vVec1(idxPeak), 'kx', 'MarkerSize', 10);
        end
        if isempty(indSpikes)
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'Location', 'SouthOutside')
        elseif isSpontaneous            % spontaneous spikes
            plot(tVec0(indSpikes(1)), vVec0(indSpikes(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike', ...
                'Location', 'SouthOutside')
        else                    % spontaneous spike or noise
            plot(tVec0(indSpikes(1)), vVec0(indSpikes(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike or noise', ...
                'Location', 'SouthOutside')
        end
    end
    plot(tVec0(idxPeakStartOld), vVec3(idxPeakStartOld), 'k*');
    plot(tVec0(idxPeakEndOld), vVec3(idxPeakEndOld), 'k*');
    plot(tVec0(idxPeakStart), vVec3(idxPeakStart), 'r*');
    plot(tVec0(idxPeakEnd), vVec3(idxPeakEnd), 'r*');
%{
    ax = gca;
    if isSpontaneous
        text((19/20)*ax.XLim(1) + (1/20)*ax.XLim(2), ...
            (1/20)*ax.YLim(1) + (19/20)*ax.YLim(2), ...
            'Spontaneous spike!');
    end
%}
    xlim([tVec0(idxPeakStart), tVec0(idxPeakEnd)]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(['Burst analysis for ', fileBaseTitle, ', original trace']);
    figname = fullfile(outFolder, outSubDirs{3}, [fileBase, '_burstanalysis', '.png']);
    saveas(h, figname);
    hold off;
    % close(h);

    % Reset default legend autoupdate
    set(groot, 'defaultLegendAutoUpdate', defaultLegendAutoUpdateOrig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fileBasesToOverride, idxMissedLtsByOrder, idxMissedLtsByShape, ...
    idxSpikesPerBurstIncorrect, idxLooksLikeMissedLts, ...
    idxLooksLikeLtsNotByProminence, ...
    idxLooksLikeLtsNotByNarrowness, ...
    idxNoiseInTrace, idxSpontLtsOrBurst, ...
    idxWideLtsCouldBeNoise] = m3ha_find_files_to_override
%% Finds traces to override

%% Hard-coded parameters
% Relative path of special cases directory under the home directory
specialCasesDir = fullfile('data_dclamp', 'take4', 'special_cases');

% Folder names to override LTS decisions
prefixToOverride = 'OVERRIDE_';
strsToOverride = {
    'Missed_LTS_by_order', ...          % LTS peak not the first below threshold
    'Missed_LTS_by_shape', ...          % LTS peak before first action potential
    'Spikes_per_burst_incorrect', ...   % Peak bounds are not correct after redefining it by maxNoise
    'Noise_in_trace', ...               % Clearly not LTSs even though detected by the algorithm
    'Spontaneous_LTSs_or_bursts', ...   % Clearly not evoked LTSs even though detected by the algorithm
    'Wide_LTS_could_be_noise', ...      % Not LTSs by vote even though detected by the algorithm
    'Looks_like_LTS_not_by_narrowness', ... % LTSs by vote even though deemed too narrow by the algorithm
    'Looks_like_LTS_not_by_prominence', ... % LTSs by vote even though deemed too small by the algorithm
    'Looks_like_missed_LTS'};           % LTSs by vote even though not detected by the algorithm
subDirsToOverride = strcat(prefixToOverride, strsToOverride);

%% Find all file base names to override
% Create full paths to directories to override
pathsToOverride = fullfile(m3ha_locate_homedir, specialCasesDir, ...
                            subDirsToOverride);

% Find all file base names to override
fileBasesToOverride = ...
    cellfun(@(x) all_file_bases('Directory', x, 'Extension', 'png'), ...
            pathsToOverride, 'UniformOutput', false);

% Remove '_scaled' from base names
fileBasesToOverride = cellfun(@(x) replace(x, '_scaled', ''), ...
                        fileBasesToOverride, 'UniformOutput', false);

% Record the indices of the directories to override
[idxMissedLtsByOrder, idxMissedLtsByShape, ...
    idxSpikesPerBurstIncorrect, idxLooksLikeMissedLts, ...
    idxLooksLikeLtsNotByProminence, ...
    idxLooksLikeLtsNotByNarrowness, ...
    idxNoiseInTrace, idxSpontLtsOrBurst, ...
    idxWideLtsCouldBeNoise] = ...
    argfun(@(x) find_in_strings(x, strsToOverride), ...
            'Missed_LTS_by_order', 'Missed_LTS_by_shape', ...
            'Spikes_per_burst_incorrect', 'Looks_like_missed_LTS', ...
            'Looks_like_LTS_not_by_prominence', ...
            'Looks_like_LTS_not_by_narrowness', ...
            'Noise_in_trace', 'Spontaneous_LTSs_or_bursts', ...
            'Wide_LTS_could_be_noise');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

idxPeak3StartOrig = zeros(nPeak3s, 1);
idxPeak3EndOrig = zeros(nPeak3s, 1);
% Detect spikes in original trace within the peak
idxPeak3StartOrig(iPeak) = idxPeak3StartThis;
idxPeak3EndOrig(iPeak) = idxPeak3EndThis;

idxPeak3StartOrig(iPeak3Sel) = idxPeak3Start(iPeak3Sel);               % beginning index of LTS peak in terms of v0
idxPeak3EndOrig(iPeak3Sel) = idxPeak3End(iPeak3Sel);                 % ending index of LTS peak in terms of v0

dvVec3 = diff(vVec3)./diff(tVec0);
d2vVec3 = diff(dvVec3Sm)./diff(tVec0(2:end));

if isempty(siMs) && ~isempty(tVec0s)
    % Compute sampling interval(s) in ms
    siMs = compute_sampling_interval(tVec0s);
elseif isempty(tVec0s) && ~isempty(siMs)
    % Create time vector(s)
    tVec0s = create_time_vectors(nSamples, 'SamplingIntervalMs', siMs, ...
                                    'TimeUnits', 'ms');
else
    error('One of siMs and tVec0s must be provided!');
end

spikeFrequency = 1000 * (numel(indSpikes)-1) / ( tVec0(indSpikes(end)) - tVec0(indSpikes(1)) );
% spike frequency (Hz) is 
% (# of spikes - 1)/(time between first and last spike)

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
