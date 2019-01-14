function varargout = parse_LTS (tVec0, vVec0, varargin)
%% Find, plot and classify the most likely low-threshold spike (LTS) candidate in a voltage trace
% Usage: [parsedParams, parsedData] = parse_LTS (tVec0, vVec0, varargin)

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
% Arguments:    
%       tVec0       - original time vector, must be a column vector in ms
%                   must be a numeric array
%       vVec0       - original voltage vector, must be a column vector in mV
%                   must be a numeric array with same length as tVec0
%       varargin    - 'StimStartMs': time of stimulation start (ms)
%                   must be a positive scalar
%                   default == TODO: first significant deflection
%                   - 'NoiseWindowMsOrMaxNoise': maximum noise in mV 
%                                           or baseline window in ms
%                                           if numel == 1, maxNoise;
%                                           if numel == 2, noiseWindowMs; 
%                   must be a numeric vector
%                   default == computed from the range [0, stimStartMs]
%                   - 'MinPeakDelayMs': minimum peak delay (ms) TODO: replace iPeakt
%                               after the end of the pulse
%                   must be a nonnegative scalar
%                   default == 0 ms
%                   - 'SearchWindowMs': window to search for LTS (ms)
%                   must be within range of tVec0
%                   default == [sStimStartMs + minPeakDelayMs, 
%                               tVec0(end) - medfiltWindowMs], where medfiltWindowMs is 30 ms
%                   - 'PlotFlag': whether to plot traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileBase': base of filename (without extension)
%                   must be a string scalar or a character vector
%                   default == 'unnamed'
%                   - 'tVec2': time vector for resampling
%                   must be a numeric array & within range of tVec0
%                   default == siMsRes*(round(tVec0(1)/siMsRes):round(tVec0(end)/siMsRes))'
%                   - 'vVec1': median-filtered voltage vector
%                   must be a numeric array with same length as vVec0
%                   default == medfilt1(vVec0, round(medfiltWindowMs/siMs))
%                   - 'vVec2': voltage vector after resampling
%                   must be a numeric array with same length as tVec2
%                   default == interp1(tVec0, vVec1, tVec2, 'linear')
%                   - 'vVec3': moving-average-filtered vVec1
%                   must be a numeric array with same length as vVec0
%                   default == smooth(vVec1, round(smoothWindowSamples))
%
% Requires:
%       cd/argfun.m
%       cd/check_subdir.m
%       cd/compute_sampling_interval.m
%       cd/all_filebases.m
%       cd/find_in_strings.m
%       cd/m3ha_locate_homedir.m
%
% Used by:

% File History:
% 2019-01-13 Adapted from find_LTS.m

%% Parameters used for data analysis
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

%% Subdirectories in outFolder for placing figures
outSubDirs = {'vtraces', 'LTSanalysis', 'burstanalysis', ...
            'vtraces_scaled', 'gray_area_traces', 'LTScouldbemissed'};

%% Traces to override
% Relative path of special cases directory under the home directory
specialCasesDir = fullfile('data_dclamp', 'take4', 'special_cases';

% Folder names to override LTS decisions
prefixToOverride = 'OVERRIDE_';
tracesToOverrideStrs = {
    'Missed_LTS_by_order', ...          % LTS peak not the first below threshold
    'Missed_LTS_by_shape', ...          % LTS peak before first action potential
    'Spikes_per_burst_incorrect', ...   % Peak bounds are not correct after redefining it by maxNoise
    'Noise_in_trace', ...               % Clearly not LTSs even though detected by the algorithm
    'Spontaneous_LTSs_or_bursts', ...   % Clearly not evoked LTSs even though detected by the algorithm
    'Wide_LTS_could_be_noise', ...      % Not LTSs by vote even though detected by the algorithm
    'Looks_like_LTS_not_by_narrowness', ... % LTSs by vote even though deemed too narrow by the algorithm
    'Looks_like_LTS_not_by_prominence', ... % LTSs by vote even though deemed too small by the algorithm
    'Looks_like_missed_LTS'};           % LTSs by vote even though not detected by the algorithm
subDirsToOverride = strcat(prefixToOverride, tracesToOverrideStrs);

%% Default values for optional arguments
stimStartMsDefault = [];        % set later
noiseWindowMsORmaxNoiseDefault = [];    % set later
minPeakDelayMsDefault = 0;
searchWindowDefault = [];       % set later
plotFlagDefault = false;
outFolderDefault = pwd;
fileBaseDefault = 'unnamed';
tVec2Default = [];              % set later
vVec1Default = [];              % set later
vVec2Default = [];              % set later
vVec3Default = [];              % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'tVec0', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec0 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'vVec0', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec0 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'NoiseWindowMsOrMaxNoise', noiseWindowMsORmaxNoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'SearchWindow', searchWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'tVec2', tVec2Default, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec2 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVec1', vVec1Default, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec1 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVec2', vVec2Default, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec2 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVec3', vVec3Default, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec3 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, tVec0, vVec0, varargin{:});
noiseWindowMsORmaxNoise = iP.Results.NoiseWindowMsOrMaxNoise;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
searchWindow = iP.Results.SearchWindow;
plotFlag = iP.Results.PlotFlag;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
tVec2 = iP.Results.tVec2;
vVec1 = iP.Results.vVec1;
vVec2 = iP.Results.vVec2;
vVec3 = iP.Results.vVec3;

% Check relationships between arguments
if ~isequal(numel(tVec0), numel(vVec0))
    error('Time and voltage vectors do not have the same length!');
end

%% Preparation
% Compute the sampling interval in ms
siMs = compute_sampling_interval(tVec0);

% Compute the base of the time vector
tBase = tVec0(1) - siMs;

% Decide whether to detect holding potential
if stimStartMs < tBase + baseWidthMs
    fprintf(['stimStartMs must be at least %g ms after tBase ', ...
                'for baseline voltage computation!\n'], baseWidthMs);
    fprintf('Actual holding potential will be empty!\n');
    actVhold = [];
    computeActVholdFlag = false;
else
    computeActVholdFlag = true;
end

% Decide on whether to compute maximum noise
if length(noiseWindowMsORmaxNoise) == 2
    noiseWindowMs = noiseWindowMsORmaxNoise;
    computeMaxNoiseFlag = true;
else
    maxNoise = noiseWindowMsORmaxNoise;
    computeMaxNoiseFlag = false;
end

% Create default search window
if isempty(searchWindow)
    searchWindow = [stimStartMs + minPeakDelayMs, tVec0(end) - medfiltWindowMs];    
end

%% Find all file base names to override
% Count the number of subdirectories to override
nDirs = numel(subDirsToOverride);

% Locate home directory
homeDirectory = m3ha_locate_homedir;

% Create full file path to special_cases directory
specialCasesPath = fullfile(homeDirectory, specialCasesDir);

% Find all file base names to override
traces_to_override = all_filebases(specialCasesPath, subDirsToOverride, 'png');

% Remove '_scaled' from base names
for k = 1:nDirs
    for t = 1:numel(traces_to_override{k})
        traces_to_override{k}{t} = replace(traces_to_override{k}{t}, '_scaled', '');
    end
end

%% Display standard output header
% fprintf('FINDING LTSs for %s ...\n', fileBase);
% fprintf('Sampling interval == %g ms\n', siMs);

% Compute things in samples
% Note: span of smooth() can't be zero
smoothWindowSamples = round(smoothWindowMs / siMs);
if smoothWindowSamples == 0
    smoothWindowSamples = 1;
end
medfiltWindowSamples = round(medfiltWindowMs/siMs);

%% Reorganize data if needed
% Median filter voltage traces to get rid of spikes
if isempty(vVec1)
    vVec1 = medfilt1(vVec0, round(medfiltWindowMs/siMs));
end

% Resample all traces at 1 kHz if needed
if isempty(tVec2)
    tVec2 = siMsRes*(round(tVec0(1)/siMsRes):round(tVec0(end)/siMsRes))';
end

% Interpolate voltage vector if needed
if isempty(vVec2)
    vVec2 = interp1(tVec0, vVec1, tVec2, 'linear');
end

% Moving-average-filter median-filtered traces for finding narrowest voltage peaks
if isempty(vVec3)
    vVec3 = smooth(vVec1, smoothWindowSamples);
end

%% Current start and peak info
idxStimStart = find(tVec0 >= stimStartMs, 1);  % Index of current application 
idxSearchStart = find(tVec0 >= stimStartMs + minPeakDelayMs, 1);   % Index of current peak
% fprintf('Time of current application == %g ms\n', stimStartMs);
% fprintf('Index of current application == %d\n', idxStimStart);
% fprintf('Time of search start == %g ms\n', stimStartMs + minPeakDelayMs);
% fprintf('Index of current peak == %d\n', idxSearchStart);

%% Find and plot LTSs
% Initialize vectors
isSpontaneous = false;      % whether it's a spontaneous spike
isOverridden = false;       % whether the algorithm is overridden
allspi = [];        % all spike indices (could be burst or spontaneous spike)
spikesPerPeak = 0;  % all peaks are initially assumed to have no spikes

% Find indices for calculating maxNoise
if computeMaxNoiseFlag
    idxNoiseBegin = find(tVec0 >= noiseWindowMs(1), 1);
    idxNoiseEnd = find(tVec0 <= noiseWindowMs(2), 1, 'last');
    indNoise = idxNoiseBegin:idxNoiseEnd;
end
if computeActVholdFlag
    indBase = idxStimStart - fliplr(1:round(baseWidthMs/siMs));    % indices for calculating baseline voltage
end
ind3_begin = find(tVec0 >= searchWindow(1), 1);
ind3_begin = max(ind3_begin, idxSearchStart);            % first index for LTS search must be after current peak
ind3_end = find(tVec0 <= searchWindow(2), 1, 'last');     % last index for LTS search

% Calculate maxNoise
if computeMaxNoiseFlag
    maxNoise = 4*std(vVec3(indNoise));    % 4*standard deviation of values of the median-filtered then 
                                        % moving-average-filtered trace is the maximum noise
                                        % Assuming a Gaussian distribution of noise, should contain 95.45%
end
% fprintf('Maximum noise == %g mV\n', maxNoise);

% Calculate baseline voltage (holding potential)
if computeActVholdFlag
    actVhold = mean(vVec1(indBase));         % baseline voltage (actual holding potential)
                                            % Previously uses vVec0
    % fprintf('Actual holding potential == %g mV\n', actVhold);
end

% Set up 2nd derivative of median-filtered voltage trace
dvvec3 = diff(vVec3)./diff(tVec0);              % differentiate vVec3, the median-filtered 
                                                %   then smoothed voltage trace
dvvec3_sm = smooth(dvvec3, smoothWindowSamples);          % smooth out dvvec3
ddvvec3 = diff(dvvec3_sm)./diff(tVec0(2:end));  % differentiate dvvec3_sm
ind3 = ind3_begin:ind3_end;
% fprintf('Finding peak in the following window: [%g %g]\n', siMs*ind3_begin, siMs*ind3_end);
v3 = vVec3(ind3);               % voltage vector of interest for detecting LTS & calculating LTS 2ndder
v0 = vVec0(ind3);               % voltage vector of interest for detecting spikes (ap)
v1 = vVec1(ind3);               % voltage vector of interest for calculating LTS amplitude, prominence, width
v_begin = ind3_begin;
ddv3 = ddvvec3(ind3 - 1);       % 2nd derivative vector of interest
                                % two differentiations: 
                                % shift right (don't -1) than shift left (-1)

% Find all voltage peaks, locate peak bounds, compute most negative 2nd derivatives, detect spikes
[vpeak_a3, vpeak_i3, vpeak_w3, vpeak_p3] = findpeaks(v3);    % find all voltage peaks
npks = length(vpeak_a3);
if npks == 0    % no local maximums exist
    [vpeak_a3, vpeak_i3] = max(v3);
    vpeak_w3 = 0;
    vpeak_p3 = 0;
    npks = 1;
end
vpeak_lb = zeros(npks, 1);
vpeak_ub = zeros(npks, 1);
vpeak_2der = zeros(npks, 1);
v0_pk = cell(npks, 1);
v0_pk_begin = zeros(npks, 1);
v0_pk_end = zeros(npks, 1);
spi = cell(npks, 1);
sp1sti = zeros(npks, 1);
for p = 1:npks
    % Find peak lower bounds
    if vpeak_i3(p) < 3            % findpeaks will not work for < 3 points
        vpeak_lb(p) = 1;
    else
        [~, ind4] = findpeaks(-flipud(v3(1:vpeak_i3(p))), ...
                'MinPeakProminence', 0, 'NPeaks', 1);
            % first minimum to the left: flip and invert, 
            % then find first voltage peak
        if isempty(ind4)
            vpeak_lb(p) = 1;
        else
            vpeak_lb(p) = (vpeak_i3(p) + 1) - ind4;
        end
    end
    % Find peak upper bounds
    if vpeak_i3(p) > length(v3) - 2    % findpeaks will not work for < 3 points
        vpeak_ub(p) = length(v3);
    else
        [~, ind5] = findpeaks(-v3(vpeak_i3(p):end), ...
                'MinPeakProminence', 0, 'NPeaks', 1);
            % first minimum to the right: invert, 
            % then find first voltage peak
        if isempty(ind5)
            vpeak_ub(p) = length(v3);
        else
            vpeak_ub(p) = (vpeak_i3(p) - 1) + ind5;
        end
    end
    % Find most negative 2nd derivative over the entire peak
    if vpeak_lb(p) == vpeak_ub(p)
        vpeak_2der(p) = ddv3(vpeak_lb(p));
    else
        vpeak_2der(p) = min(ddv3(vpeak_lb(p):vpeak_ub(p)));
    end
    % Detect spikes in original trace within the peak
    v0_pk{p} = v0(vpeak_lb(p):vpeak_ub(p));        % take only the peak part
    v0_pk_begin(p) = vpeak_lb(p);
    v0_pk_end(p) = vpeak_ub(p);
    [pspikes_a, pspikes_i] = findpeaks(v0_pk{p});    % find all "spikes" within the peak
    stemp1 = find(pspikes_a > spThr);        % action potentials must be greater than threshold
    if ~isempty(stemp1)
        % Record spike indices relative to v0 or v3
        spi{p} = (v0_pk_begin(p) - 1) + pspikes_i(stemp1);
        % Record index of first spike relative to v0 or v3
        sp1sti(p) = spi{p}(1);
    end
end

% Algorithm for detecting whether there is an LTS and classifying other peaks:
% (1) Eliminate all voltage peaks with prominence < maxNoise
%     then find the first voltage peak with second derivative below threshold.
%     If it has spikes, 
%        it is an LTS only if the first spike occurs before the peak.
%
%    Old: Eliminate all voltage peaks with prominence < 1 mV 
%     or with the first spike occurring after the peak, 
%     then find the first voltage peak with second derivative below threshold. 
%
%     Old: Eliminate all voltage peaks with prominence < 1 mV 
%     or with the first spike occurring after the peak, 
%    then find the narrowest voltage peak; It's an LTS only 
%    if the second derivative is below threshold.
% (2) If doesn't exist, 
%    eliminate all voltage peaks with prominence < maxNoise, 
%     then find the narrowest voltage peak.
% (3) If still doesn't exist, find the narrowest voltage peak.

% Find all peaks with prominence greater than maximum noise
ptemp1 = find(vpeak_p3 > maxNoise);            % peak #s with prom > maxNoise         
                            % (previously 1 mV, now different for each trace)

% Find all peaks with the second derivative reaching (<=) threshold
ptemp2 = find(vpeak_2der <= ltsThr);

% Find all peaks that are LTS candidates by threshold
ptemp3 = intersect(ptemp1, ptemp2);            % peak #s that are LTSs by threshold

if isempty(ptemp1)        % Condition (3)
    % Not LTS by prominence
    ltsPeakValue = NaN;
    ltsPeakTime = NaN;
    maxSlopeTime = NaN;
    maxSlopeValue = NaN;

    % Find the narrowest voltage peak
    [peak2ndDer, psel3] = min(vpeak_2der);        % find the narrowest voltage peak
    peakClass = 1;
    peakClassLabel = ['not LTS: peak prominence ', ...
                num2str(vpeak_p3(psel3)), ...
                ' mV <= maximum noise ', ...
                num2str(maxNoise), ' mV'];
    
    % Override the algorithm if 3/4 experts think it's an LTS
    ind_LlLnbp = find_in_strings('Looks_like_LTS_not_by_prominence', tracesToOverrideStrs);
    if ismember(fileBase, traces_to_override{ind_LlLnbp})
        isOverridden = true;
        peakClass = 5;
        peakClassLabel = ['LTS by overrule: peak prominence ', ...
                num2str(vpeak_p3(psel3)), ...
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
    [peak2ndDer, pp1] = min(vpeak_2der(ptemp1));
    psel3 = ptemp1(pp1);
    peakClass = 2;
    peakClassLabel = ['not LTS: 2nd derivative ', ...
                num2str(peak2ndDer), ...
                ' V^2/s^2 > LTS threshold ', ...
                num2str(ltsThr), ' V^2/s^2'];

    % Override the algorithm if 3/4 experts think it's an LTS
    ind_LlLnbn = find_in_strings('Looks_like_LTS_not_by_narrowness', tracesToOverrideStrs);
    if ismember(fileBase, traces_to_override{ind_LlLnbn})
        isOverridden = true;
        peakClass = 5;
        peakClassLabel = ['LTS by overrule: 2nd derivative ', ...
                    num2str(peak2ndDer), ...
                    ' V^2/s^2 > LTS threshold ', ...
                    num2str(ltsThr), ' V^2/s^2'];
    end

    % There were a few cases (see 'Looks_like_missed_LTS') where the correct LTS was missed by prominence %%% CHECK FOR MORE
    ind_LlmL = find_in_strings('Looks_like_missed_LTS', tracesToOverrideStrs);
    if ismember(fileBase, traces_to_override{ind_LlmL})
        psel3 = ptemp2(1);        % the first narrow peak regardless of prominence is the correct one in this case
        peak2ndDer = vpeak_2der(psel3);    % update peak2ndDer
        isOverridden = true;
        peakClass = 5;
        peakClassLabel = ['LTS by overrule: peak prominence ', ...
                num2str(vpeak_p3(psel3)), ...
                ' mV <= maximum noise ', ...
                num2str(maxNoise), ' mV'];
    end
else                % Condition (1)
    % Select the first peak that is an LTS candidate by prominence & 2nd der
    % There is one case (see 'Missed_LTS_by_order') where the second peak is the correct LTS
    ind_MLbo = find_in_strings('Missed_LTS_by_order', tracesToOverrideStrs);
    if ismember(fileBase, traces_to_override{ind_MLbo})
        psel3 = ptemp3(2);        % for F092710_0006_25
    else
        psel3 = ptemp3(1);
    end
    peak2ndDer = vpeak_2der(psel3);        % update peak2ndDer

    % Check whether it's a spontaneous spike
    %     (Based on following observation of shape:
    %         LTS:         first spike occurs before "LTS" peak on mfmaf trace, 
    %            except in four cases (see 'Missed_LTS_by_shape'),
    %             where the first spike occurred after the peak
    %     spontaneous spike:     first spike occurs after "LTS" peak on mfmaf trace)
    sp2pk_i = round(minSp2PkTime/siMs);    % obsolete: currently set to 0
    ind_MLbs = find_in_strings('Missed_LTS_by_shape', tracesToOverrideStrs);
    if ~ismember(fileBase, traces_to_override{ind_MLbs}) ...
        && sp1sti(psel3) ~= 0 ...
        && vpeak_i3(psel3) < sp1sti(psel3) + sp2pk_i
        % Not LTS by shape
        ltsPeakValue = NaN;
        ltsPeakTime = NaN;
        maxSlopeTime = NaN;
        maxSlopeValue = NaN;

        isSpontaneous = false;
        peakClass = 3;
        peakClassLabel = ['not LTS: peak index ', ...
                    num2str(vpeak_i3(psel3)), ...
                    ' < index of first spike ', ...
                    num2str(sp1sti(psel3))];
    end
end
% fprintf('Selected peak 2nd derivative == %g V^2/s^2\n', peak2ndDer);

% Record width and prominence of selected peak
[vpeak_a1, vpeak_i1, vpeak_w1, vpeak_p1] = findpeaks(v1);    % find all voltage peaks for v1
npks1 = length(vpeak_a1);
if npks1 == 0    % no local maximums exist
    [vpeak_a1, vpeak_i1] = max(v1);
    vpeak_w1 = 0;
    vpeak_p1 = 0;
    npks = 1;
end
[~, psel1] = min(abs(vpeak_i1 - vpeak_i3(psel3)));        % find peak of v1 that corresponds to selected peak of v3
peakProm = vpeak_p1(psel1);            % this is the "minimum vertical distance that the signal must 
                        % descend on either side of the peak before either climbing back 
                        % to a level higher than the peak or reaching an endpoint" in mV
peakWidth = vpeak_w1(psel1) * siMs;        % the width of the peak at half-prominence in ms

% fprintf('Selected peak prominence == %g mV\n', peakProm);
% fprintf('Selected peak width == %g ms\n', peakWidth);

% Record indices relative to tVec0 or vVec0
npi = (v_begin - 1) + vpeak_i1(psel1);        % narrowest peak index, relative to median-filtered voltage trace
peakTime = (npi - idxStimStart) * siMs;        % narrowest peak time (delay), relative to median-filtered voltage trace
% fprintf('Selected peak time == %g ms\n', peakTime);

% The following may be changed later for bursts
np_lbi = (v_begin - 1) + vpeak_lb(psel3);    % narrowest peak lower bound index, relative to moving-average-filtered voltage trace
np_ubi = (v_begin - 1) + vpeak_ub(psel3);    % narrowest peak upper bound index, relative to moving-average-filtered voltage trace

% Save old peak boundary indices for plotting
np_lbi_old = np_lbi;
np_ubi_old = np_ubi;

% Record spike indices relative to vVec0 & count spikes per peak
if ~isempty(spi{psel3})
    allspi = (v_begin - 1) + spi{psel3};    % spike indices
    spikesPerPeak = length(spi{psel3});     % count spikes per peak
end

% Find LTS amplitude, delay, maximum slope delay and value
if (~isempty(ptemp3) && ~isSpontaneous) ...    % either peakClass hasn't been classified yet
    || isOverridden                         % or it's an 'LTS by overrule'
    ltsPeakValue = vpeak_a1(psel1);           % LTS amplitude, use median-filtered voltage trace
    ltsPeakTime = peakTime;                 % LTS delay, use median-filtered voltage trace
    sp_thr = ltsPeakValue + spThrRelLts;   % threshold for spike detection (mV)

    % Find approximate max slope based off of v3
    spacing_size = round(slopeSpacing/siMs);                           % number of indices apart to calculate slope
    v3pk_left = v3(vpeak_lb(psel3):vpeak_ub(psel3) - spacing_size);     % the left points for slope calculation
    v3pk_right = v3(vpeak_lb(psel3) + spacing_size:vpeak_ub(psel3));    % the right points for slope calculation
    v3pk_slopes = (v3pk_right - v3pk_left)/slopeSpacing;               % all slopes in V/s
    [maxslopeval_appr, ~] = max(v3pk_slopes);                           % find approximate maximum slope
%    [maxSlopeValue, temp_ind] = max(v3pk_slopes);                        % find maximum slope

%{
    % Find max slope based off of v1
    spacing_size = round(slopeSpacing/siMs);                           % number of indices apart to calculate slope
    v1pk_left = v1(vpeak_lb(psel3):vpeak_ub(psel3) - spacing_size);     % the left points for slope calculation
    v1pk_right = v1(vpeak_lb(psel3) + spacing_size:vpeak_ub(psel3));    % the right points for slope calculation
    v1pk_slopes = (v1pk_right - v1pk_left)/slopeSpacing;               % all slopes in V/s
    [maxSlopeValue, temp_ind] = max(v1pk_slopes);                         % find maximum slope
%}

    % Moving-average-filter median-filtered traces for calculating maximum slope
    mafw3 = mafw3Dv/maxslopeval_appr;      % width in ms for the moving average filter for finding slopes
    ndp_mafw3 = round(mafw3/siMs);
    if ndp_mafw3 == 0                       % span of smooth() can't be zero
        ndp_mafw3 = 1;
    end
    v4 = smooth(v1, ndp_mafw3);             % voltage vector of interest for calculating maxslope

    % Find max slope based off of v4
    spacing_size = round(slopeSpacing/siMs);                       % number of indices apart to calculate slope
    v4pk_left = v4(vpeak_lb(psel3):vpeak_ub(psel3) - spacing_size); % the left points for slope calculation
    v4pk_right = v4(vpeak_lb(psel3) + spacing_size:vpeak_ub(psel3));% the right points for slope calculation
    v4pk_slopes = (v4pk_right - v4pk_left)/slopeSpacing;           % all slopes in V/s
    [maxSlopeValue, temp_ind] = max(v4pk_slopes);                     % find maximum slope

    maxslopeind = temp_ind + spacing_size + (np_lbi - 1);           % the maxslope index is the right point
    maxSlopeTime = (maxslopeind - idxStimStart) * siMs;               % delay in ms of maximum slope after IPSC starts
    peakfeature_label = ['max slope = ', num2str(maxSlopeValue), ' V/s; ', ...
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
    if isempty(allspi)          % no spikes detected, not a "burst"
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

        ind_Spbi = find_in_strings('Spikes_per_burst_incorrect', tracesToOverrideStrs);
        if ~ismember(fileBase, traces_to_override{ind_Spbi}) ...
            % Re-detect spikes by redefining peak bounds using maxNoise as MinPeakProminence
            % Update peak lower bounds
            if vpeak_i3(psel3) < 3                  % findpeaks will not work for < 3 points
                vpeak_lb(psel3) = 1;
            else
                [~, ind4] = findpeaks(-flipud(v3(1:vpeak_i3(psel3))), ...
                        'MinPeakProminence', maxNoise, 'NPeaks', 1);
                    % first minimum to the left: flip and invert, 
                    % then find first voltage peak with prominence 
                    % >= maxNoise (previously 2 mV, now different for each trace)
                if isempty(ind4)
                    vpeak_lb(psel3) = 1;
                else
                    vpeak_lb(psel3) = (vpeak_i3(psel3) + 1) - ind4;
                end
            end
            % Update peak upper bounds
            if vpeak_i3(psel3) > length(v3) - 2    % findpeaks will not work for < 3 points
                vpeak_ub(psel3) = length(v3);
            else
                [~, ind5] = findpeaks(-v3(vpeak_i3(psel3):end), ...
                        'MinPeakProminence', maxNoise, 'NPeaks', 1);
                    % first minimum to the right: invert, 
                    % then find first voltage peak with prominence 
                    % >= maxNoise (previously 2 mV, now different for each trace)
                if isempty(ind5)
                    vpeak_ub(psel3) = length(v3);
                else
                    vpeak_ub(psel3) = (vpeak_i3(psel3) - 1) + ind5;
                end
            end
            % Detect spikes in original trace within the peak
            v0_pk{psel3} = v0(vpeak_lb(psel3):vpeak_ub(psel3)); % take only the peak part
            v0_pk_begin(psel3) = vpeak_lb(psel3);               % beginning index of LTS peak in terms of v0
            v0_pk_end(psel3) = vpeak_ub(psel3);                 % ending index of LTS peak in terms of v0
            [pspikes_a, pspikes_i] = findpeaks(v0_pk{psel3});   % find all "spikes" within the peak
            stemp1 = find(pspikes_a > sp_thr);          % action potentials must be greater than actual threshold
                                % Note: this is a different threshold than before,
                                %    so could potentially change the classification of the peak
            if ~isempty(stemp1)
                % Record spike indices relative to v0 or v3
                spi{psel3} = (v0_pk_begin(psel3) - 1) + pspikes_i(stemp1);
            else
                spi{psel3} = [];
            end

            % Update peak bound indices relative to tVec0 or vVec0
            np_lbi = (v_begin - 1) + vpeak_lb(psel3);       % narrowest peak lower bound index
            np_ubi = (v_begin - 1) + vpeak_ub(psel3);       % narrowest peak upper bound index
        end

        % Two cases: spikes are still found or not
        if ~isempty(spi{psel3})
            % Update spike indices and spikes per peak
            allspi = (v_begin - 1) + spi{psel3};            % spike indices
            spikesPerPeak = length(allspi);                 % count spikes per peak

            % Record time of first and last spike
            firstSpikeTime = (allspi(1) - idxStimStart) * siMs;
            lastSpikeTime = (allspi(end) - idxStimStart) * siMs;
            spiketime_label = ['spikes; first = ', num2str(firstSpikeTime), ...
                                ' ms; last = ', num2str(lastSpikeTime), ' ms'];

            % Spikes per burst is the the same spikes per peak but not zero
            spikesPerBurst = spikesPerPeak;                 % spikes per burst

            % Define the burst region
            vvec0_b = v0_pk{psel3};                         % voltage trace of burst
            bspk1i = sp1sti(psel3) - v0_pk_begin(psel3) + 1;% index of first spike in vvec0_b

            % Find first minimum to the left of first spike
            if bspk1i < 3                       % findpeaks will not work for < 3 points
                spk_lb = 1;
                amp6 = -vvec0_b(1);
            else
                [amp6, ind6] = findpeaks(-flipud(vvec0_b(1:bspk1i)), 'NPeaks', 1);
                if isempty(ind6)
                    spk_lb = 1;
                    amp6 = -vvec0_b(1);
                else
                    spk_lb = (bspk1i + 1) - ind6;
                end
            end
            % Find first minimum to the right of first spike
            if bspk1i > length(vvec0_b) - 2     % findpeaks will not work for < 3 points
                spk_ub = length(vvec0_b);
                amp7 = -vvec0_b(end);
            else
                [amp7, ind7] = findpeaks(-vvec0_b(bspk1i:end), 'NPeaks', 1);
                if isempty(ind7)
                    spk_ub = length(vvec0_b);
                    amp7 = -vvec0_b(end);
                else
                    spk_ub = (bspk1i - 1) + ind7;
                end
            end
            amp8 = max([-amp6 -amp7]);          % take the higher of the two minimums

            % Find the last index lower than the base of the first spike in terms of vvec0_b
            ind8 = find(vvec0_b(1:bspk1i) < amp8, 1, 'last');

            % Find the burst onset time (delay)
            if isempty(ind8) || ind8 < 4
                % Burst onset index is the beginning of the LTS peak in terms of vVec0
                bon_i = (v_begin - 1) + v0_pk_begin(psel3);
            else                        % trace to differentiate must have at least 4 elements
                % Find the corresponding median-filtered then smoothed trace for the burst region
                vvec4_b = v4(v0_pk_begin(psel3):v0_pk_end(psel3));
                
                % Find the index (in terms of vvec4_b) of the maximum of the 3rd derivative of
                %   the voltage trace in between the start of burst and the last point before the first spike
                %   Add 3 indices to account for the loss of an element on the left after each diff()
                [~, ind9] = max(diff(diff(diff(vvec4_b(1:ind8)))) / (siMs^3));
                max3rdderi = ind9 + 3;
                
                % Burst onset index is in terms of vVec0
                bon_i = (v_begin - 1) + (v0_pk_begin(psel3) - 1) + max3rdderi;
            end
            burstTime = (bon_i - idxStimStart) * siMs;    % burst onset time (delay)
            bursttime_label = ['burst onset; delay = ', num2str(burstTime), ' ms'];

            % The computed spike threshold is the voltage at burst onset
            spikeThreshold = vVec0(bon_i);

            % Find the maximum spike amplitude, minimum spike amplitude, 
            %   spike frequency and spike adaptation
            if length(allspi) >= 1
                maxSpikeAmp = max(vVec0(allspi));
                minSpikeAmp = min(vVec0(allspi));
                if length(allspi) >= 2
                    spikeFrequency = 1000 * (length(allspi)-1) / ( tVec0(allspi(end)) - tVec0(allspi(1)) );
                                    % spike frequency (Hz) is 
                                    % (# of spikes - 1)/(time between first and last spike)
                    if length(allspi) >= 3
                        spikeAdaptation = 100 * (tVec0(allspi(end)) - tVec0(allspi(end-1))) ...
                                            / (tVec0(allspi(2)) - tVec0(allspi(1)));
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
            allspi = [];
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
[indNoiseInTrace, indSpontLtsOrBurst, indWideLtsCouldBeNoise] = ...
    argfun(@(x) find_in_strings(x, tracesToOverrideStrs), ...
            'Noise_in_trace', 'Spontaneous_LTSs_or_bursts', ...
            'Wide_LTS_could_be_noise');

if ismember(fileBase, traces_to_override{indNoiseInTrace}) ...
    || ismember(fileBase, traces_to_override{indSpontLtsOrBurst}) ...
    || ismember(fileBase, traces_to_override{indWideLtsCouldBeNoise})
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

%% Plot figures
if plotFlag
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
        leftslopepoint_y = min(vVec3(np_ubi), vVec3(maxslopeind) - slopeSegYHalf);
        rightslopepoint_y = min(vVec3(npi), vVec3(maxslopeind) + slopeSegYHalf);
        leftslopepoint_x = (leftslopepoint_y - vVec3(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);
        rightslopepoint_x = (rightslopepoint_y - vVec3(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);
%}

        leftslopepoint_y = min(vVec1(np_ubi), vVec1(maxslopeind) - slopeSegYHalf);
        rightslopepoint_y = min(vVec1(npi), vVec1(maxslopeind) + slopeSegYHalf);
        leftslopepoint_x = (leftslopepoint_y - vVec1(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);
        rightslopepoint_x = (rightslopepoint_y - vVec1(maxslopeind)) / maxSlopeValue + tVec0(maxslopeind);

        % Compute info for peak prominence line segment
        pkprom_x = tVec0(npi);            % time value (ms) of peak prominence
        pkprom_y1 = vVec1(npi) - peakProm;    % voltage value (mV) of bottom of peak
        pkprom_y2 = vVec1(npi);            % voltage value (mV) of top of peak

        % Compute info for peak width line segment
        halfprom_y = vVec1(npi)- peakProm/2;            % voltage value (mV) at half prominence
        [~, temp_ind] = min(abs(vVec1(np_lbi:npi) - halfprom_y));    % Find the index on the rising phase 
                                        % of the peak whose voltage value 
                                        % is closest to half prominence
        pkwidth_leftind = temp_ind + np_lbi - 1;        % peak width segment left index relative to vVec1
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
            plot(tVec0(npi), vVec1(npi), 'go', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'ro', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                % line for maxslope
        plot(tVec0(bon_i), vVec0(bon_i), 'g>', 'MarkerSize', 10);       % triangle for burst onset time
        plot(tVec0(allspi), vVec0(allspi), 'gx', 'MarkerSize', 10);     % crosses for spikes
        legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakfeature_label, bursttime_label, spiketime_label, ...
            'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
        line([tVec0(1), tVec0(bon_i)], [spikeThreshold, spikeThreshold], ...
                'LineStyle', ':');                                      % line for spike threshold
        line([tVec0(1), tVec0(end)], [maxSpikeAmp, maxSpikeAmp], 'LineStyle', '--');    % line for maxSpikeAmp
        line([tVec0(1), tVec0(end)], [minSpikeAmp, minSpikeAmp], 'LineStyle', '--');    % line for minSpikeAmp
        text(.05, .15, ['Spike Frequency: ', num2str(spikeFrequency), ' Hz'], 'Units', 'normalized');
        text(.05, .1, ['Spike Adaptation: ',  num2str(spikeAdaptation), '%'], 'Units', 'normalized');
        text(.05, .05, ['Spike Threshold: ',  num2str(spikeThreshold), 'mV'], 'Units', 'normalized');
    elseif ~isnan(ltsPeakTime)            % LTS without bursts
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), vVec1(npi), 'bo', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(npi), vVec1(npi), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                             % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                              % line for maxslope
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
        legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakfeature_label, 'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
    else                        % not LTS
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(npi), vVec1(npi), 'kx', 'MarkerSize', 10);
        end
        if isempty(allspi)            % noise
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'Location', 'SouthOutside')
        elseif isSpontaneous            % spontaneous spikes
            plot(tVec0(allspi(1)), vVec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike', ...
                'Location', 'SouthOutside')
        else                    % spontaneous spike or noise
            plot(tVec0(allspi(1)), vVec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike or noise', ...
                'Location', 'SouthOutside')
        end
    end
    if computeActVholdFlag
        plot(tVec0(idxBaseStart), vVec1(idxBaseStart), 'g>');
        plot(tVec0(idxBaseEnd), vVec1(idxBaseEnd), 'y<');
    end
    plot(tVec0(np_lbi_old), vVec3(np_lbi_old), 'k*');
    plot(tVec0(np_ubi_old), vVec3(np_ubi_old), 'k*');
    plot(tVec0(np_lbi), vVec3(np_lbi), 'r*');
    plot(tVec0(np_ubi), vVec3(np_ubi), 'r*');
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
    xlimits = searchWindow;
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
    xlimits = searchWindow;
    h = figure(5001);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'LTS analysis, moving-average-filtered trace');
    clf(h);
    subplot(3,1,1) % voltage trace
    plot(tVec0, vVec3, 'r-', 'LineWidth', 0.5); hold on
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), vVec1(npi), 'go', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'gd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(npi), vVec1(npi), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'md', 'MarkerSize', 8);
        end
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), vVec1(npi), 'bo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'bd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(npi), vVec1(npi), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), vVec1(maxslopeind), 'md', 'MarkerSize', 8);
        end
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(npi), vVec1(npi), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tVec0(np_lbi), vVec3(np_lbi), 'r*');
    plot(tVec0(np_ubi), vVec3(np_ubi), 'r*');
    plot(tVec0(np_lbi_old), vVec3(np_lbi_old), 'k*');
    plot(tVec0(np_ubi_old), vVec3(np_ubi_old), 'k*');
    xlim(xlimits);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(['LTS analysis for ', fileBaseTitle, ', moving-average-filtered trace']);
    subplot(3,1,2) % 1st derivative of voltage trace
    plot(tVec0(2:end), dvvec3, 'k-', 'LineWidth', 0.5); hold on
    plot(tVec0(2:end), dvvec3_sm, 'r-', 'LineWidth', 0.5);
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), dvvec3(npi), 'go', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvvec3(maxslopeind), 'gd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), dvvec3(npi), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvvec3(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(npi), dvvec3(npi), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvvec3(maxslopeind), 'md', 'MarkerSize', 8);
        end
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), dvvec3(npi), 'bo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvvec3(maxslopeind), 'bd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), dvvec3(npi), 'ro', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvvec3(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(npi), dvvec3(npi), 'mo', 'MarkerSize', 10);
            plot(tVec0(maxslopeind), dvvec3(maxslopeind), 'md', 'MarkerSize', 8);
        end
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), dvvec3(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(npi), dvvec3(npi), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tVec0(np_lbi), dvvec3(np_lbi), 'r*');
    plot(tVec0(np_ubi), dvvec3(np_ubi), 'r*');
    plot(tVec0(np_lbi_old), dvvec3(np_lbi_old), 'k*');
    plot(tVec0(np_ubi_old), dvvec3(np_ubi_old), 'k*');
    xlim(xlimits);
    xlabel('Time (ms)');
    ylabel('dV/dT');
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
    legend('unsmoothed', 'smoothed');
    subplot(3,1,3) % 2nd derivative of voltage trace
    plot(tVec0(3:end), ddvvec3, 'k-', 'LineWidth', 0.5); hold on
    line(xlimits, [ltsThr ltsThr], 'Color', 'r', ...
        'LineStyle', '--', 'LineWidth', 0.5);                       % mark LTS threshold
    line(xlimits, [ltsThrAlt ltsThrAlt], 'Color', 'b', ...
        'LineStyle', '--', 'LineWidth', 0.5);                       % mark alt LTS threshold
    xlim(xlimits);
    ax = gca;
    ylimits = get(ax, 'YLim');
    line([tVec0(ind3_begin) tVec0(ind3_begin)], ...
        ylimits, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);% mark start of search window
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), ddvvec3(npi), 'go', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr        % in "gray area"
            plot(tVec0(npi), ddvvec3(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr         % LTS by overrule
            plot(tVec0(npi), ddvvec3(npi), 'mo', 'MarkerSize', 10);
        end
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), ddvvec3(npi), 'bo', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr        % in "gray area"
            plot(tVec0(npi), ddvvec3(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr         % LTS by overrule
            plot(tVec0(npi), ddvvec3(npi), 'mo', 'MarkerSize', 10);
        end
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr        % in "gray area"
            plot(tVec0(npi), ddvvec3(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(npi), ddvvec3(npi), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tVec0(np_lbi), ddvvec3(np_lbi), 'r*');
    plot(tVec0(np_ubi), ddvvec3(np_ubi), 'r*');
    plot(tVec0(np_lbi_old), ddvvec3(np_lbi_old), 'k*');
    plot(tVec0(np_ubi_old), ddvvec3(np_ubi_old), 'k*');
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
    plot(tVec0(np_lbi:np_ubi), vVec0(np_lbi:np_ubi), 'b-', 'LineWidth', 0.5); hold on
    plot(tVec2, vVec2, 'g.', 'LineWidth', 0.5);
    plot(tVec0(np_lbi:np_ubi), vVec3(np_lbi:np_ubi), 'r-', 'LineWidth', 0.5);
    if ~isnan(ltsPeakTime) && ~isnan(burstTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), vVec1(npi), 'go', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(npi), vVec1(npi), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                             % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                              % line for maxslope
        plot(tVec0(bon_i), vVec0(bon_i), 'g>', 'MarkerSize', 10);       % triangle for burst onset time
        plot(tVec0(allspi), vVec0(allspi), 'gx', 'MarkerSize', 10);     % crosses for spikes
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
        legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakfeature_label, bursttime_label, spiketime_label, ...
            'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
        line([tVec0(1), tVec0(bon_i)], [spikeThreshold, spikeThreshold], ...
                'LineStyle', ':');                                      % line for spike threshold
        line([tVec0(1), tVec0(end)], [maxSpikeAmp, maxSpikeAmp], 'LineStyle', '--');    % line for maxSpikeAmp
        line([tVec0(1), tVec0(end)], [minSpikeAmp, minSpikeAmp], 'LineStyle', '--');    % line for minSpikeAmp
        text(.05, .7, ['Spike Frequency: ', num2str(spikeFrequency), ' Hz'], 'Units', 'normalized');
        text(.05, .65, ['Spike Adaptation: ',  num2str(spikeAdaptation), '%'], 'Units', 'normalized');
        text(.05, .6, ['Spike Threshold: ',  num2str(spikeThreshold), 'mV'], 'Units', 'normalized');
    elseif ~isnan(ltsPeakTime)
        if peak2ndDer <= ltsThrAlt
            plot(tVec0(npi), vVec1(npi), 'bo', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndDer > ltsThr        % LTS by overrule
            plot(tVec0(npi), vVec1(npi), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                % line for maxslope
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
        legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
            peakClassLabel, peakfeature_label, 'Location', 'SouthOutside')
%        plot(tVec0(maxslopeind), vVec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tVec0(maxslopeind), vVec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakProm
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakWidth
    else
        if peak2ndDer > ltsThrAlt ...
            && peak2ndDer <= ltsThr    % in "gray area"
            plot(tVec0(npi), vVec1(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tVec0(npi), vVec1(npi), 'kx', 'MarkerSize', 10);
        end
        if isempty(allspi)
    %% TODO: Fix the legend: Label the plots and use 'Displayname' and legend(subset, ...)
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'Location', 'SouthOutside')
        elseif isSpontaneous            % spontaneous spikes
            plot(tVec0(allspi(1)), vVec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike', ...
                'Location', 'SouthOutside')
        else                    % spontaneous spike or noise
            plot(tVec0(allspi(1)), vVec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakClassLabel, 'spontaneous spike or noise', ...
                'Location', 'SouthOutside')
        end
    end
    plot(tVec0(np_lbi_old), vVec3(np_lbi_old), 'k*');
    plot(tVec0(np_ubi_old), vVec3(np_ubi_old), 'k*');
    plot(tVec0(np_lbi), vVec3(np_lbi), 'r*');
    plot(tVec0(np_ubi), vVec3(np_ubi), 'r*');
%{
    ax = gca;
    if isSpontaneous
        text((19/20)*ax.XLim(1) + (1/20)*ax.XLim(2), ...
            (1/20)*ax.YLim(1) + (19/20)*ax.YLim(2), ...
            'Spontaneous spike!');
    end
%}
    xlim([tVec0(np_lbi), tVec0(np_ubi)]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(['Burst analysis for ', fileBaseTitle, ', original trace']);
    figname = fullfile(outFolder, outSubDirs{3}, [fileBase, '_burstanalysis', '.png']);
    saveas(h, figname);
    hold off;
    % close(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
