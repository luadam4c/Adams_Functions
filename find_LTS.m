function [actVhold, maxnoise, peaktime, peak2ndder, peakprom, peakwidth, peakclass, spikesperpeak, ltspeaktime, ltspeakval, maxslopetime, maxslopeval, bursttime, spikesperburst, spikethreshold, firstspiketime, lastspiketime, maxspikeamp, minspikeamp, spikefrequency, spikeadaptation, LTScouldbemissed] = find_LTS (tvec0, vvec0, istart, ipeakt, hrangeORmaxnoise, ltswin, plotflag, outfolder, filebase, tvec2, vvec1, vvec2, vvec3)
%% Find, plot and classify the most likely low-threshold spike (LTS) candidate in a voltage trace
% Usage: [actVhold, maxnoise, peaktime, peak2ndder, peakprom, peakwidth, peakclass, spikesperpeak, ltspeaktime, ltspeakval, maxslopetime, maxslopeval, bursttime, spikesperburst, spikethreshold, firstspiketime, lastspiketime, maxspikeamp, minspikeamp, spikefrequency, spikeadaptation, LTScouldbemissed] = find_LTS (tvec0, vvec0, istart, ipeakt, hrangeORmaxnoise, ltswin, plotflag, outfolder, filebase, tvec2, vvec1, vvec2, vvec3)
% Arguments:    
%       tvec0       - original time vector, must be a column vector in ms
%       vvec0       - original voltage vector, must be a column vector in mV
%                   must be a numeric array with same length as tvec0
%       istart      - the time point (in ms) from which peaktime, ltspeaktime, bursttime are calculated, e.g., 1000
%                   usually the time of current application
%                   must be within range of tvec0 and >= tvec0(1) + baseline width (20 ms)
%       ipeakt      - the time point (in ms) after which the LTS is detected, e.g., 1300
%                   usually the time of current peak
%                   must be within range of tvec0 and >= istart
%       hrangeORmaxnoise - if length == 2, hrange; if length == 1, maxnoise
%       hrange      - Range (in ms) in which maximum noise is calculated, e.g. [200 1000]
%                   usually a time interval before istart
%                   must be within range of tvec0
%       maxnoise    - maximum noise in mV, e.g., 1
%                   must be >= 0
%       ltswin      - (opt) Window (in ms) in which the LTS would lie, e.g. [1000 7960]
%                   must be within range of tvec0
%                   default == [ipeakt tvec0(end)-mfw3], where mfw3 is 30 ms
%       plotflag    - (opt) whether to plot traces or not
%                   must be 0 or 1
%                   default == 0
%       outfolder   - (opt) directory to place outputs, e.g. '/media/adamX/m3ha/data_dclamp/take4/test_sweeps/'
%                   must be a directory
%                   default == pwd
%       filebase    - (opt) base of filename (without extension), e.g. 'A100110_0008_18'
%                   must be a char array
%                   default == 'unnamed'
%       tvec2       - (opt) time vector for resampling
%                   must be a numeric array & within range of tvec0
%                   default == rsims*(round(tvec0(1)/rsims):round(tvec0(end)/rsims))'
%       vvec1       - (opt) median-filtered voltage vector
%                   must be a numeric array with same length as vvec0
%                   default == medfilt1(vvec0, round(mfw3/sims))
%       vvec2       - (opt) voltage vector after resampling
%                   must be a numeric array with same length as tvec2
%                   default == interp1(tvec0, vvec1, tvec2, 'linear')
%       vvec3       - (opt) moving-average-filtered vvec1
%                   must be a numeric array with same length as vvec0
%                   default == smooth(vvec1, round(ndp_mafw2))
%
% Requires:    
%        cd/check_subdir.m
%        cd/find_filebases.m
%        cd/find_ind_str_in_cell.m
%
% Used by:
%        /home/Matlab/Adams_Functions/find_LTSs_many_sweeps
%        /media/adamX/m3ha/data_dclamp/test_sweep.m
%        /media/adamX/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m
%        /media/shareX/share/Adam/Sample_files_from_Katie/test_sweeps.m
%
%
% 2016-09-22 Adapted from dclampDataExtractor.m
% 2016-10-05 Accounted for the case that tvec0 doesn't start from 0
% 2016-10-05 Made maxnoise an argument that could replace hrange
% 2016-10-13 Now plots maxslope with a diamond
% 2016-10-13 Changed maxnoise so that it's 2*standard deviation of the baseline data instead of the maximal range
% 2016-10-14 Changed spike threshold from -30 mV to -43 mV
% 2016-10-14 Made tvec2, vvec1, vvec2, vvec3 optional arguments for efficiency
% 2016-10-16 Added back close(h)
% 2016-10-17 Changed spike threshold to 5 mV above LTS peak (Added sp_thr_rel_lts = 5)
% 2016-10-17 Changed maxnoise so that it's 4*standard deviation of the baseline data instead of 2*standard deviation
% 2016-10-17 Commented out close(h) and put close all in for loop instead
% 2016-10-18 Changed spike threshold again to 15 mV above LTS peak (sp_thr_rel_lts = 15)
% 2016-10-18 Changed sp_thr_init to -45 mV from -43 mV 
% 2016-10-19 Changed spike threshold again to 10 mV above LTS peak (sp_thr_rel_lts = 10) 
% 2016-10-19 Now prints "spontaneous spike or noise" for peaks of class 1 & 2 with spikes
% 2016-10-19 Fixed issue with not updating detected spikes (see D092710_0006_5) spi{psel3} = []
% 2016-10-19 max slope now displays a value when plotted
% 2016-10-19 Changed the y axis of vtraces_scaled and plotted median-filtered trace with thicker line 
% 2016-10-24 BT - Added lines for minimum and maximum spike amplitudes
% 2016-10-27 BT - Changed spike frequency to be computed from peak to peak
% 2016-10-27 Replace each directory with directories{k}
% 2016-10-31 BT - Added line for maxslope
% 2016-11-01 AL - Fixed peakwidth so that it's in ms
% 2016-11-01 BT - Added lines for peakprom and peakwidth
% 2016-11-02 AL - Moved check directories to check_subdir.m
% 2016-11-02 AL - Changed maxslope_label to peakfeature_label
% 2016-11-17 BT & AL - Changed peaktime, peakprom, peakwidth to be computed 
%               from the median-filtered trace instead of the moving-average-filtered trace
% 2016-11-17 BT - Changed the plots of peaktime, peakprom, peakwidth accordingly
% 2016-11-21 BT - Changed maxslope to be computed from the median-filtered trace 
%               instead of the moving-average-filtered trace
% 2016-11-29 BT - Added spacing_size; Changed maxslope to be computed over 1 ms instead of over 0.1 ms
% 2016-11-29 AL - Added slope_spacing; Changed maxslope to be computed with vector manipulation for performance
% 2016-11-29 AL - Made sure maxslope, peakprom, peakwidth was plotted on LTSs without bursts as well
% 2016-11-29 AL - Changed the limits of the maxslope line segment so it would show up at the right place
% 2016-11-29 AL - Added mafw3, vvec4 & v4 for finding maximum slope
% 2016-11-30 AL - Added slopesegyhalf for plotting maximum slope
% 2016-11-30 AL - Added mafw3_dv, changed mafw3, v4 to depend on the maximum slope found from v3
% 2017-01-16 AL - Accounted for the condition that npks == 0 (no local maximums exist)
% 2017-02-15 AL - Added traces_to_override
% 2017-02-16 AL - Fixed 'Missed_LTS_by_order', 'Missed_LTS_by_shape'
% 2017-02-16 AL - Changed peakclass to 1~9 (from 1~7) and peakclass_label to nine items
% 2017-02-16 AL - Fixed 'Noise_in_trace'
% 2017-03-20 AL - Fixed 'Spikes_per_burst_incorrect', 'Spontaneous_LTSs_or_bursts', 'Wide_LTS_could_be_noise'
% 2017-03-20 AL - Fixed 'Looks_like_LTS_not_by_narrowness', 'Looks_like_LTS_not_by_prominence'
% 2017-03-22 BT - Added LTScouldbemissed
% 2017-03-22 AL - Fixed 'Looks_like_missed_LTS' partially
% 2017-05-19 AL - Made sure ndp_mafw2 & ndp_mafw3 are always positive
% 2017-07-28 AL - Fixed xlimits of burstanalysis plots
% 2017-12-21 AL - Changed tabs to spaces
% 2017-12-21 AL - Added spike threshold, first spike time, last spike time
% 2017-12-21 AL - Changed the definition of burst onset time
% 2017-12-21 AL - Added bursttime_label and spiketime_label
%

%% Parameters used for data analysis
mfw3 = 30;          % width in ms for the median filter for spikes (voltage traces)
mafw2 = 30;         % width in ms for the moving average filter for finding narrowest voltage peaks
blw = 20;           % width in ms for calculating baseline voltage (holding potential)
lts_thr = -0.0023;  % 2nd derivative in V^2/s^2 below which defines an LTS peak
lts_thr_alt = -0.0081823;   % 2nd derivative in V^2/s^2 above which is the "gray area"
sp_thr_init = -45;  % Initial amplitude threshold in mV for detecting a spike 
                    % Will be changed later when LTS peak amplitude is found
                    % 2016-10-14 the highest LTS peak without bursts is -43.8867 mV
                    % 2016-10-18 Note: the highest LTS peak without bursts is actually -48.0006 mV
                    %           the previous value was for a spontaneous LTS, but because of median-filtering,
                    %           -45 mV is probably a safer threshold
sp_thr_rel_lts = 10;% Relative amplitude threshold in mV for detecting a spike above an LTS 
                    % 2016-10-19 The smallest relative amplitude of an action potential riding above the LTS
                    %           is probably between 10~11 mV (see B091010_0006_19)
rsims = 1;          % resampling interval in ms (1 kHz)
sp2pk_t = 0;        % minimum time from the first spike to the peak of the LTS (ms)
slope_spacing = 1;  % spacing in ms used to calculate slope
mafw3_dv = 3;       % voltage change in mV corresponding to the moving average filter window for finding slopes
%mafw3 = 5;         % width in ms for the moving average filter for finding slopes
slopesegyhalf = 5;  % how much voltage difference (mV) to plot maxslope line segment below maxslope point

%% Directories for placing figures
directories = {'/vtraces/', '/LTSanalysis/', '/burstanalysis/', '/vtraces_scaled/', '/gray_area_traces/', '/LTScouldbemissed/'};

%% Traces to override
dir_special_cases = '/data_dclamp/take4/special_cases/';        % under homedirectory
traces_to_override_strs = {
    'Missed_LTS_by_order', ...          % LTS peak not the first below threshold
    'Missed_LTS_by_shape', ...          % LTS peak before first action potential
    'Spikes_per_burst_incorrect', ...   % Peak bounds are not correct after redefining it by maxnoise
    'Noise_in_trace', ...               % Clearly not LTSs even though detected by the algorithm
    'Spontaneous_LTSs_or_bursts', ...   % Clearly not evoked LTSs even though detected by the algorithm
    'Wide_LTS_could_be_noise', ...      % Not LTSs by vote even though detected by the algorithm
    'Looks_like_LTS_not_by_narrowness', ... % LTSs by vote even though deemed too narrow by the algorithm
    'Looks_like_LTS_not_by_prominence', ... % LTSs by vote even though deemed too small by the algorithm
    'Looks_like_missed_LTS'};           % LTSs by vote even though not detected by the algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 5
    error('Not enough input arguments, type ''help find_LTS'' for usage');
elseif isempty(tvec0) || isempty(vvec0) || isempty(istart) || isempty(ipeakt) || isempty(hrangeORmaxnoise) 
    error('First five inputs cannot be empty!');
elseif ~isnumeric(tvec0) || ~isnumeric(vvec0) || ~isnumeric(istart) || ~isnumeric(ipeakt) || ~isnumeric(hrangeORmaxnoise) 
    error('First five inputs must be numbers or numeric arrays!');
elseif size(tvec0, 2) > 1 || size(vvec0, 2) > 1 
    error('tvec0 and vvec0 must be column vectors!');
elseif ~isequal(size(tvec0, 1), size(vvec0, 1))
    error('Time and voltage vectors do not have the same length!');
end

%% Extract info from data
sims = tvec0(2) - tvec0(1);            % sampling interval in ms
tbase = tvec0(1) - sims;

%% Check more arguments
if istart < tbase + blw || istart > tvec0(end) || ipeakt < istart || ipeakt > tvec0(end)
    error('istart or ipeakt out of range!');
elseif length(hrangeORmaxnoise) > 2
    error('hrange or maxnoise has too many components!');
elseif length(hrangeORmaxnoise) == 2 && (hrangeORmaxnoise(1) < tbase || hrangeORmaxnoise(2) > tvec0(end))
    error('hrange out of range!');
elseif length(hrangeORmaxnoise) == 1 && hrangeORmaxnoise < 0
    error('maxnoise out of range!');
elseif nargin >= 6 && length(ltswin) < 2
    error('ltswin must have a start and an end!');
elseif nargin >= 6 && (ltswin(1) < tbase || ltswin(2) > tvec0(end))
    error('ltswin out of range!');
elseif nargin >= 7 && (plotflag ~= 1 && plotflag ~= 0 && plotflag ~= false && plotflag ~= true)
    error('plotflag out of range!');
elseif nargin >= 8 && ~ischar(outfolder)
    error('outfolder must be a character array!');
elseif nargin >= 9 && ~ischar(filebase)
    error('filebase must be a character array!');
elseif nargin >= 10 && ~isnumeric(tvec2)
    error('tvec2 must be a numeric array!');
elseif nargin >= 10 && (tvec2(1) < tbase || tvec2(end) > tvec0(end) + sims)
    error('tvec2 out of range!');
elseif nargin >= 11 && ~isnumeric(vvec1)
    error('vvec1 must be a numeric array!');
elseif nargin >= 11 && ~isequal(size(vvec0, 1), size(vvec1, 1))
    error('vvec1 must have the same length as vvec0!');
elseif nargin >= 12 && ~isnumeric(vvec2)
    error('vvec2 must be a numeric array!');
elseif nargin >= 12 && ~isequal(size(tvec2, 1), size(vvec2, 1))
    error('Time and voltage vectors do not have the same length!');
elseif nargin >= 13 && ~isnumeric(vvec3)
    error('vvec3 must be a numeric array!');
elseif nargin >= 13 && ~isequal(size(vvec0, 1), size(vvec3, 1))
    error('vvec3 must have the same length as vvec0!');
end

%% Set defaults for optional arguments
if length(hrangeORmaxnoise) == 2
    hrange = hrangeORmaxnoise;
else
    maxnoise = hrangeORmaxnoise;
end
if nargin < 6
    ltswin = [ipeakt tvec0(end)-mfw3];    
end
if nargin < 7
    plotflag = 0;
end
if nargin < 8
    outfolder = pwd;
end
if nargin < 9
    filebase = 'unnamed';
end

%% Locate home directory
if exist('/media/adamX/m3ha/', 'dir') == 7
    homedirectory = '/media/adamX/m3ha/';
elseif exist('/scratch/al4ng/m3ha/', 'dir') == 7
    homedirectory = '/scratch/al4ng/m3ha/';
else
    error('Valid homedirectory does not exist!');
end

%% Find all filebases to override
ndirs = numel(traces_to_override_strs);         % number of subdirectories in the special_cases directory
traces_to_override_subdirs = cell(1, ndirs);    % store the names of the subdirectories
for k = 1:ndirs
    traces_to_override_subdirs{k} = ['/OVERRIDE_', traces_to_override_strs{k}, '/'];
end
dir_special_cases_full = fullfile(homedirectory, dir_special_cases);    % full file path to special_cases directory
traces_to_override = find_filebases(dir_special_cases_full, traces_to_override_subdirs, 'png');

% Remove '_scaled' from base names
for k = 1:ndirs
    for t = 1:numel(traces_to_override{k})
        traces_to_override{k}{t} = strrep(traces_to_override{k}{t}, '_scaled', '');
    end
end

%% Display standard output header
% fprintf('FINDING LTSs for %s ...\n', filebase);
% fprintf('Sampling interval == %g ms\n', sims);

%% Reorganize data if needed
% Median filter voltage traces to get rid of spikes
if nargin < 10
    vvec1 = medfilt1(vvec0, round(mfw3/sims));
end
% Resample all traces at 1 kHz for fitting use (Christine's traces)
if nargin < 11
    tvec2 = rsims*(round(tvec0(1)/rsims):round(tvec0(end)/rsims))';     % resampled time vector
end
if nargin < 12
    vvec2 = interp1(tvec0, vvec1, tvec2, 'linear');                     % interpolated voltage vector
end
% Moving-average-filter median-filtered traces for finding narrowest voltage peaks
ndp_mafw2 = round(mafw2/sims);
if ndp_mafw2 == 0                       % span of smooth() can't be zero
    ndp_mafw2 = 1;
end
if nargin < 13
    vvec3 = smooth(vvec1, ndp_mafw2);
end

%% Current start and peak info
istart_ind = find(tvec0 >= istart, 1);  % Index of current application 
ipeak_ind = find(tvec0 >= ipeakt, 1);   % Index of current peak
% fprintf('Time of current application == %g ms\n', istart);
% fprintf('Index of current application == %d\n', istart_ind);
% fprintf('Time of current peak == %g ms\n', ipeakt);
% fprintf('Index of current peak == %d\n', ipeak_ind);

%% Find and plot LTSs
% Initialize vectors
isspontaneous = 0;                      % Whether it's a spontaneous spike
isoverridden = 0;                       % Whether the algorithm is overridden
LTScouldbemissed = 0;                   % Whether the LTS might have been missed
allspi = [];                            % All spike indices (could be burst or spontaneous spike)
spikesperpeak = 0;                      % All peaks are initially assumed to have no spikes

% Find indices
if length(hrangeORmaxnoise) == 2
    hr_ind_begin = find(tvec0 >= hrange(1), 1);
    hr_ind_end = find(tvec0 <= hrange(2), 1, 'last');
    hr_ind = hr_ind_begin:hr_ind_end;   % indices for calculating maxnoise
end
bl_ind = istart_ind - fliplr(1:round(blw/sims));    % indices for calculating baseline voltage
ind3_begin = find(tvec0 >= ltswin(1), 1);
ind3_begin = max(ind3_begin, ipeak_ind);            % first index for LTS search must be after current peak
ind3_end = find(tvec0 <= ltswin(2), 1, 'last');     % last index for LTS search

% Calculate maxnoise
if length(hrangeORmaxnoise) == 2
    maxnoise = 4*std(vvec3(hr_ind));    % 4*standard deviation of values of the median-filtered then 
                                        % moving-average-filtered trace is the maximum noise
                                        % Assuming a Gaussian distribution of noise, should contain 95.45%
end
% fprintf('Maximum noise == %g mV\n', maxnoise);

% Calculate baseline voltage (holding potential)
actVhold = mean(vvec1(bl_ind));         % baseline voltage (actual holding potential)
                                        % Previously uses vvec0
% fprintf('Actual holding potential == %g mV\n', actVhold);

% Set up 2nd derivative of median-filtered voltage trace
dvvec3 = diff(vvec3)./diff(tvec0);              % differentiate vvec3, the median-filtered 
                                                %   then smoothed voltage trace
dvvec3_sm = smooth(dvvec3, ndp_mafw2);          % smooth out dvvec3
ddvvec3 = diff(dvvec3_sm)./diff(tvec0(2:end));  % differentiate dvvec3_sm
ind3 = ind3_begin:ind3_end;
% fprintf('Finding peak in the following window: [%g %g]\n', sims*ind3_begin, sims*ind3_end);
v3 = vvec3(ind3);               % voltage vector of interest for detecting LTS & calculating LTS 2ndder
v0 = vvec0(ind3);               % voltage vector of interest for detecting spikes (ap)
v1 = vvec1(ind3);               % voltage vector of interest for calculating LTS amplitude, prominence, width
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
    stemp1 = find(pspikes_a > sp_thr_init);        % action potentials must be greater than threshold
    if ~isempty(stemp1)
        % Record spike indices relative to v0 or v3
        spi{p} = (v0_pk_begin(p) - 1) + pspikes_i(stemp1);
        % Record index of first spike relative to v0 or v3
        sp1sti(p) = spi{p}(1);
    end
end

% Algorithm for detecting whether there is an LTS and classifying other peaks:
% (1) Eliminate all voltage peaks with prominence < maxnoise
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
%    eliminate all voltage peaks with prominence < maxnoise, 
%     then find the narrowest voltage peak.
% (3) If still doesn't exist, find the narrowest voltage peak.

% Find all peaks with prominence greater than maximum noise
ptemp1 = find(vpeak_p3 > maxnoise);            % peak #s with prom > maxnoise         
                            % (previously 1 mV, now different for each trace)

% Find all peaks with the second derivative reaching (<=) threshold
ptemp2 = find(vpeak_2der <= lts_thr);

% Find all peaks that are LTS candidates by threshold
ptemp3 = intersect(ptemp1, ptemp2);            % peak #s that are LTSs by threshold

if isempty(ptemp1)        % Condition (3)
    % Not LTS by prominence
    ltspeakval = NaN;
    ltspeaktime = NaN;
    maxslopetime = NaN;
    maxslopeval = NaN;

    % Find the narrowest voltage peak
    [peak2ndder, psel3] = min(vpeak_2der);        % find the narrowest voltage peak
    peakclass = 1;
    peakclass_label = ['not LTS: peak prominence ', ...
                num2str(vpeak_p3(psel3)), ...
                ' mV <= maximum noise ', ...
                num2str(maxnoise), ' mV'];
    
    % Override the algorithm if 3/4 experts think it's an LTS
    ind_LlLnbp = find_ind_str_in_cell('Looks_like_LTS_not_by_prominence', traces_to_override_strs);
    if ismember(filebase, traces_to_override{ind_LlLnbp})
        isoverridden = 1;
        peakclass = 5;
        peakclass_label = ['LTS by overrule: peak prominence ', ...
                num2str(vpeak_p3(psel3)), ...
                ' mV <= maximum noise ', ...
                num2str(maxnoise), ' mV'];
    end
elseif isempty(ptemp3) ...    % Condition (2)
    % Not LTS by narrowness
    ltspeakval = NaN;
    ltspeaktime = NaN;
    maxslopetime = NaN;
    maxslopeval = NaN;

    % Find the narrowest voltage peak with prominence greater than maximum noise
    [peak2ndder, pp1] = min(vpeak_2der(ptemp1));
    psel3 = ptemp1(pp1);
    peakclass = 2;
    peakclass_label = ['not LTS: 2nd derivative ', ...
                num2str(peak2ndder), ...
                ' V^2/s^2 > LTS threshold ', ...
                num2str(lts_thr), ' V^2/s^2'];

    % Override the algorithm if 3/4 experts think it's an LTS
    ind_LlLnbn = find_ind_str_in_cell('Looks_like_LTS_not_by_narrowness', traces_to_override_strs);
    if ismember(filebase, traces_to_override{ind_LlLnbn})
        isoverridden = 1;
        peakclass = 5;
        peakclass_label = ['LTS by overrule: 2nd derivative ', ...
                    num2str(peak2ndder), ...
                    ' V^2/s^2 > LTS threshold ', ...
                    num2str(lts_thr), ' V^2/s^2'];
    end

    % There were a few cases (see 'Looks_like_missed_LTS') where the correct LTS was missed by prominence %%% CHECK FOR MORE
    ind_LlmL = find_ind_str_in_cell('Looks_like_missed_LTS', traces_to_override_strs);
    if ismember(filebase, traces_to_override{ind_LlmL})
        psel3 = ptemp2(1);        % the first narrow peak regardless of prominence is the correct one in this case
        peak2ndder = vpeak_2der(psel3);    % update peak2ndder
        isoverridden = 1;
        peakclass = 5;
        peakclass_label = ['LTS by overrule: peak prominence ', ...
                num2str(vpeak_p3(psel3)), ...
                ' mV <= maximum noise ', ...
                num2str(maxnoise), ' mV'];
    end
else                % Condition (1)
    % Select the first peak that is an LTS candidate by prominence & 2nd der
    % There is one case (see 'Missed_LTS_by_order') where the second peak is the correct LTS
    ind_MLbo = find_ind_str_in_cell('Missed_LTS_by_order', traces_to_override_strs);
    if ismember(filebase, traces_to_override{ind_MLbo})
        psel3 = ptemp3(2);        % for F092710_0006_25
    else
        psel3 = ptemp3(1);
    end
    peak2ndder = vpeak_2der(psel3);        % update peak2ndder

    % Check whether it's a spontaneous spike
    %     (Based on following observation of shape:
    %         LTS:         first spike occurs before "LTS" peak on mfmaf trace, 
    %            except in four cases (see 'Missed_LTS_by_shape'),
    %             where the first spike occurred after the peak
    %     spontaneous spike:     first spike occurs after "LTS" peak on mfmaf trace)
    sp2pk_i = round(sp2pk_t/sims);    % obsolete: currently set to 0
    ind_MLbs = find_ind_str_in_cell('Missed_LTS_by_shape', traces_to_override_strs);
    if ~ismember(filebase, traces_to_override{ind_MLbs}) ...
        && sp1sti(psel3) ~= 0 ...
        && vpeak_i3(psel3) < sp1sti(psel3) + sp2pk_i
        % Not LTS by shape
        ltspeakval = NaN;
        ltspeaktime = NaN;
        maxslopetime = NaN;
        maxslopeval = NaN;

        isspontaneous = 1;
        peakclass = 3;
        peakclass_label = ['not LTS: peak index ', ...
                    num2str(vpeak_i3(psel3)), ...
                    ' < index of first spike ', ...
                    num2str(sp1sti(psel3))];
    end
end
% fprintf('Selected peak 2nd derivative == %g V^2/s^2\n', peak2ndder);

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
peakprom = vpeak_p1(psel1);            % this is the "minimum vertical distance that the signal must 
                        % descend on either side of the peak before either climbing back 
                        % to a level higher than the peak or reaching an endpoint" in mV
peakwidth = vpeak_w1(psel1) * sims;        % the width of the peak at half-prominence in ms

% fprintf('Selected peak prominence == %g mV\n', peakprom);
% fprintf('Selected peak width == %g ms\n', peakwidth);

% Record indices relative to tvec0 or vvec0
npi = (v_begin - 1) + vpeak_i1(psel1);        % narrowest peak index, relative to median-filtered voltage trace
peaktime = (npi - istart_ind) * sims;        % narrowest peak time (delay), relative to median-filtered voltage trace
% fprintf('Selected peak time == %g ms\n', peaktime);

% The following may be changed later for bursts
np_lbi = (v_begin - 1) + vpeak_lb(psel3);    % narrowest peak lower bound index, relative to moving-average-filtered voltage trace
np_ubi = (v_begin - 1) + vpeak_ub(psel3);    % narrowest peak upper bound index, relative to moving-average-filtered voltage trace

% Save old peak boundary indices for plotting
np_lbi_old = np_lbi;
np_ubi_old = np_ubi;

% Record spike indices relative to vvec0 & count spikes per peak
if ~isempty(spi{psel3})
    allspi = (v_begin - 1) + spi{psel3};    % spike indices
    spikesperpeak = length(spi{psel3});     % count spikes per peak
end

% Find LTS amplitude, delay, maximum slope delay and value
if (~isempty(ptemp3) && isspontaneous ~= 1) ...    % either peakclass hasn't been classified yet
    || isoverridden                         % or it's an 'LTS by overrule'
    ltspeakval = vpeak_a1(psel1);           % LTS amplitude, use median-filtered voltage trace
    ltspeaktime = peaktime;                 % LTS delay, use median-filtered voltage trace
    sp_thr = ltspeakval + sp_thr_rel_lts;   % threshold for spike detection (mV)

    % Find approximate max slope based off of v3
    spacing_size = round(slope_spacing/sims);                           % number of indices apart to calculate slope
    v3pk_left = v3(vpeak_lb(psel3):vpeak_ub(psel3) - spacing_size);     % the left points for slope calculation
    v3pk_right = v3(vpeak_lb(psel3) + spacing_size:vpeak_ub(psel3));    % the right points for slope calculation
    v3pk_slopes = (v3pk_right - v3pk_left)/slope_spacing;               % all slopes in V/s
    [maxslopeval_appr, ~] = max(v3pk_slopes);                           % find approximate maximum slope
%    [maxslopeval, temp_ind] = max(v3pk_slopes);                        % find maximum slope

%{
    % Find max slope based off of v1
    spacing_size = round(slope_spacing/sims);                           % number of indices apart to calculate slope
    v1pk_left = v1(vpeak_lb(psel3):vpeak_ub(psel3) - spacing_size);     % the left points for slope calculation
    v1pk_right = v1(vpeak_lb(psel3) + spacing_size:vpeak_ub(psel3));    % the right points for slope calculation
    v1pk_slopes = (v1pk_right - v1pk_left)/slope_spacing;               % all slopes in V/s
    [maxslopeval, temp_ind] = max(v1pk_slopes);                         % find maximum slope
%}

    % Moving-average-filter median-filtered traces for calculating maximum slope
    mafw3 = mafw3_dv/maxslopeval_appr;      % width in ms for the moving average filter for finding slopes
    ndp_mafw3 = round(mafw3/sims);
    if ndp_mafw3 == 0                       % span of smooth() can't be zero
        ndp_mafw3 = 1;
    end
    v4 = smooth(v1, ndp_mafw3);             % voltage vector of interest for calculating maxslope

    % Find max slope based off of v4
    spacing_size = round(slope_spacing/sims);                       % number of indices apart to calculate slope
    v4pk_left = v4(vpeak_lb(psel3):vpeak_ub(psel3) - spacing_size); % the left points for slope calculation
    v4pk_right = v4(vpeak_lb(psel3) + spacing_size:vpeak_ub(psel3));% the right points for slope calculation
    v4pk_slopes = (v4pk_right - v4pk_left)/slope_spacing;           % all slopes in V/s
    [maxslopeval, temp_ind] = max(v4pk_slopes);                     % find maximum slope

    maxslopeind = temp_ind + spacing_size + (np_lbi - 1);           % the maxslope index is the right point
    maxslopetime = (maxslopeind - istart_ind) * sims;               % delay in ms of maximum slope after IPSC starts
    peakfeature_label = ['max slope = ', num2str(maxslopeval), ' V/s; ', ...
                'prominence = ', num2str(peakprom), ' mV; ', ...
                'width = ', num2str(peakwidth), ' ms'];
end

% Determine whether it's a "burst" and count action potentials (spikes)
if isnan(ltspeaktime)                % must be an "LTS" to begin with
    bursttime = NaN;
    spikesperburst = NaN;
    spikethreshold = NaN;
    firstspiketime = NaN;
    lastspiketime = NaN;
    maxspikeamp = NaN;
    minspikeamp = NaN;
    spikefrequency = NaN;
    spikeadaptation = NaN;
else
    if isempty(allspi)            % no spikes detected, not a "burst"
        bursttime = NaN;
        spikesperburst = NaN;
        spikethreshold = NaN;
        firstspiketime = NaN;
        lastspiketime = NaN;
        maxspikeamp = NaN;
        minspikeamp = NaN;
        spikefrequency = NaN;
        spikeadaptation = NaN;
        if peak2ndder <= lts_thr_alt && ~isoverridden
            peakclass = 8;
            peakclass_label = ['LTS with no burst; ', ...
                        'definite: 2nd der ', ...
                        num2str(peak2ndder), ...
                        ' V^2/s^2'];
        elseif ~isoverridden
            peakclass = 6;
            peakclass_label = ['LTS with no burst; ', ...
                        'contentious: 2nd der ', ...
                        num2str(peak2ndder), ...
                        ' V^2/s^2 > ', ...
                        'alt LTS thr ', ...
                        num2str(lts_thr_alt), ...
                        ' V^2/s^2'];
        end
    else
        if peak2ndder <= lts_thr_alt && ~isoverridden
            peakclass = 9;
            peakclass_label = ['LTS with burst; definite: ', ...
                        '2nd der ', ...
                        num2str(peak2ndder), ...
                        ' V^2/s^2'];
        elseif ~isoverridden
            peakclass = 7;
            peakclass_label = ['LTS with burst; contentious: ', ...
                        '2nd der ', ...
                        num2str(peak2ndder), ...
                        ' V^2/s^2 > ', ...
                        'alt LTS thr ', ...
                        num2str(lts_thr_alt), ...
                        ' V^2/s^2'];
        end

        ind_Spbi = find_ind_str_in_cell('Spikes_per_burst_incorrect', traces_to_override_strs);
        if ~ismember(filebase, traces_to_override{ind_Spbi}) ...
            % Re-detect spikes by redefining peak bounds using maxnoise as MinPeakProminence
            % Update peak lower bounds
            if vpeak_i3(psel3) < 3                  % findpeaks will not work for < 3 points
                vpeak_lb(psel3) = 1;
            else
                [~, ind4] = findpeaks(-flipud(v3(1:vpeak_i3(psel3))), ...
                        'MinPeakProminence', maxnoise, 'NPeaks', 1);
                    % first minimum to the left: flip and invert, 
                    % then find first voltage peak with prominence 
                    % >= maxnoise (previously 2 mV, now different for each trace)
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
                        'MinPeakProminence', maxnoise, 'NPeaks', 1);
                    % first minimum to the right: invert, 
                    % then find first voltage peak with prominence 
                    % >= maxnoise (previously 2 mV, now different for each trace)
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

            % Update peak bound indices relative to tvec0 or vvec0
            np_lbi = (v_begin - 1) + vpeak_lb(psel3);       % narrowest peak lower bound index
            np_ubi = (v_begin - 1) + vpeak_ub(psel3);       % narrowest peak upper bound index
        end

        % Two cases: spikes are still found or not
        if ~isempty(spi{psel3})
            % Update spike indices and spikes per peak
            allspi = (v_begin - 1) + spi{psel3};            % spike indices
            spikesperpeak = length(allspi);                 % count spikes per peak

            % Record time of first and last spike
            firstspiketime = (allspi(1) - istart_ind) * sims;
            lastspiketime = (allspi(end) - istart_ind) * sims;
            spiketime_label = ['spikes; first = ', num2str(firstspiketime), ...
                                ' ms; last = ', num2str(lastspiketime), ' ms'];

            % Spikes per burst is the the same spikes per peak but not zero
            spikesperburst = spikesperpeak;                 % spikes per burst

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
                % Burst onset index is the beginning of the LTS peak in terms of vvec0
                bon_i = (v_begin - 1) + v0_pk_begin(psel3);
            else                        % trace to differentiate must have at least 4 elements
                % Find the corresponding median-filtered then smoothed trace for the burst region
                vvec4_b = v4(v0_pk_begin(psel3):v0_pk_end(psel3));
                
                % Find the index (in terms of vvec4_b) of the maximum of the 3rd derivative of
                %   the voltage trace in between the start of burst and the last point before the first spike
                %   Add 3 indices to account for the loss of an element on the left after each diff()
                [~, ind9] = max(diff(diff(diff(vvec4_b(1:ind8)))) / (sims^3));
                max3rdderi = ind9 + 3;
                
                % Burst onset index is in terms of vvec0
                bon_i = (v_begin - 1) + (v0_pk_begin(psel3) - 1) + max3rdderi;
            end
            bursttime = (bon_i - istart_ind) * sims;    % burst onset time (delay)
            bursttime_label = ['burst onset; delay = ', num2str(bursttime), ' ms'];

            % The computed spike threshold is the voltage at burst onset
            spikethreshold = vvec0(bon_i);

            % Find the maximum spike amplitude, minimum spike amplitude, 
            %   spike frequency and spike adaptation
            if length(allspi) >= 1
                maxspikeamp = max(vvec0(allspi));
                minspikeamp = min(vvec0(allspi));
                if length(allspi) >= 2
                    spikefrequency = 1000 * (length(allspi)-1) / ( tvec0(allspi(end)) - tvec0(allspi(1)) );
                                    % spike frequency (Hz) is 
                                    % (# of spikes - 1)/(time between first and last spike)
                    if length(allspi) >= 3
                        spikeadaptation = 100 * (tvec0(allspi(end)) - tvec0(allspi(end-1))) ...
                                            / (tvec0(allspi(2)) - tvec0(allspi(1)));
                                    % spike adaptation (%) is last ISI/first ISI
                    else
                        spikeadaptation = NaN;
                    end
                else
                    spikefrequency = NaN;
                    spikeadaptation = NaN;
                end
            else
                maxspikeamp = NaN;
                minspikeamp = NaN;
                spikefrequency = NaN;
                spikeadaptation = NaN;
            end
        else                        % no longer a burst
            allspi = [];
            spikesperpeak = 0;
            bursttime = NaN;
            spikesperburst = NaN;
            spikethreshold = NaN;
            firstspiketime = NaN;
            lastspiketime = NaN;
            maxspikeamp = NaN;
            minspikeamp = NaN;
            spikefrequency = NaN;
            spikeadaptation = NaN;
            if peak2ndder <= lts_thr_alt && ~isoverridden
                peakclass = 8;
                peakclass_label = ['LTS with no burst; ', ...
                            'definite: 2nd der ', ...
                            num2str(peak2ndder), ...
                            ' V^2/s^2'];
            elseif ~isoverridden
                peakclass = 6;
                peakclass_label = ['LTS with no burst; ', ...
                            'contentious: 2nd der ', ...
                            num2str(peak2ndder), ...
                            ' V^2/s^2 > ', ...
                            'alt LTS thr ', ...
                            num2str(lts_thr_alt), ...
                            ' V^2/s^2'];
            end
        end
    end
end

%% Deal with false positives
ind_Nit = find_ind_str_in_cell('Noise_in_trace', traces_to_override_strs);
ind_SLob = find_ind_str_in_cell('Spontaneous_LTSs_or_bursts', traces_to_override_strs);
ind_WLcbn = find_ind_str_in_cell('Wide_LTS_could_be_noise', traces_to_override_strs);

if ismember(filebase, traces_to_override{ind_Nit}) ...
    || ismember(filebase, traces_to_override{ind_SLob}) ...
    || ismember(filebase, traces_to_override{ind_WLcbn})
    % Not LTS by overrule
    ltspeakval = NaN;
    ltspeaktime = NaN;
    maxslopetime = NaN;
    maxslopeval = NaN;
    bursttime = NaN;
    spikesperburst = NaN;
    maxspikeamp = NaN;
    minspikeamp = NaN;
    spikefrequency = NaN;
    spikeadaptation = NaN;
    isspontaneous = 0;

    % Label info according to the previously detected peakclass
    if peakclass == 6 || peakclass == 7
        peakclass_label = ['not LTS by overrule: 2nd der ', ...
                    num2str(peak2ndder), ...
                    ' V^2/s^2 < ', ...
                    'LTS threshold ', ...
                    num2str(lts_thr), ...
                    ' V^2/s^2'];
    elseif peakclass == 8 || peakclass == 9
        peakclass_label = ['not LTS by overrule: 2nd der ', ...
                    num2str(peak2ndder), ...
                    ' V^2/s^2 < ', ...
                    'alt LTS thr ', ...
                    num2str(lts_thr_alt), ...
                    ' V^2/s^2'];
    end
    
    % Make this peakclass == 4
    peakclass = 4;            
end

% fprintf('Selected peak class == %d\n', peakclass);
% fprintf('Spikes per peak == %d\n', spikesperpeak);
% fprintf('LTS peak time == %g ms\n', ltspeaktime);
% fprintf('LTS amplitude (absolute) == %g mV\n', ltspeakval);
% fprintf('LTS maximum slope time == %g ms\n', maxslopetime);
% fprintf('LTS maximum slope amplitude == %g V/s\n', maxslopeval);
% fprintf('Burst onset time == %g ms\n', bursttime);
% fprintf('Spikes per burst == %d\n\n', spikesperburst);

if plotflag

    % Check if needed directories exist in outfolder
    check_subdir(outfolder, directories);

    % For convenience
    bl_start = bl_ind(1);
    bl_end = bl_ind(end);

    % Compute info needed for plotting peak properties
    if ~isnan(ltspeaktime)
        % Compute maxslope line segment limits
%{
        leftslopepoint_y = min(vvec3(np_ubi), vvec3(maxslopeind) - slopesegyhalf);
        rightslopepoint_y = min(vvec3(npi), vvec3(maxslopeind) + slopesegyhalf);
        leftslopepoint_x = (leftslopepoint_y - vvec3(maxslopeind)) / maxslopeval + tvec0(maxslopeind);
        rightslopepoint_x = (rightslopepoint_y - vvec3(maxslopeind)) / maxslopeval + tvec0(maxslopeind);
%}

        leftslopepoint_y = min(vvec1(np_ubi), vvec1(maxslopeind) - slopesegyhalf);
        rightslopepoint_y = min(vvec1(npi), vvec1(maxslopeind) + slopesegyhalf);
        leftslopepoint_x = (leftslopepoint_y - vvec1(maxslopeind)) / maxslopeval + tvec0(maxslopeind);
        rightslopepoint_x = (rightslopepoint_y - vvec1(maxslopeind)) / maxslopeval + tvec0(maxslopeind);

        % Compute info for peak prominence line segment
        pkprom_x = tvec0(npi);            % time value (ms) of peak prominence
        pkprom_y1 = vvec1(npi) - peakprom;    % voltage value (mV) of bottom of peak
        pkprom_y2 = vvec1(npi);            % voltage value (mV) of top of peak

        % Compute info for peak width line segment
        halfprom_y = vvec1(npi)- peakprom/2;            % voltage value (mV) at half prominence
        [~, temp_ind] = min(abs(vvec1(np_lbi:npi) - halfprom_y));    % Find the index on the rising phase 
                                        % of the peak whose voltage value 
                                        % is closest to half prominence
        pkwidth_leftind = temp_ind + np_lbi - 1;        % peak width segment left index relative to vvec1
        pkwidth_x1 = tvec0(pkwidth_leftind);            % time value (ms) of left end
        pkwidth_x2 = tvec0(pkwidth_leftind) + peakwidth;    % time value (ms) of right end
    end

    % Plot voltage traces
    % fprintf('Plotting voltage trace ...\n');
    h = figure(5000);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Voltage trace');
    clf(h);
    plot(tvec0, vvec0, 'b-', 'LineWidth', 0.5); hold on;
    plot(tvec0, vvec1, 'g-', 'LineWidth', 0.5);
    plot(tvec0, vvec3, 'r-', 'LineWidth', 0.5);
    if ~isnan(ltspeaktime) && ~isnan(bursttime)    % LTS with bursts
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), vvec1(npi), 'go', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'ro', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                % line for maxslope
        plot(tvec0(bon_i), vvec0(bon_i), 'g>', 'MarkerSize', 10);       % triangle for burst onset time
        plot(tvec0(allspi), vvec0(allspi), 'gx', 'MarkerSize', 10);     % crosses for spikes
        legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
            peakclass_label, peakfeature_label, bursttime_label, spiketime_label, ...
            'Location', 'SouthOutside')
%        plot(tvec0(maxslopeind), vvec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tvec0(maxslopeind), vvec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakprom
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakwidth
        line([tvec0(1), tvec0(bon_i)], [spikethreshold, spikethreshold], ...
                'LineStyle', ':');                                      % line for spike threshold
        line([tvec0(1), tvec0(end)], [maxspikeamp, maxspikeamp], 'LineStyle', '--');    % line for maxspikeamp
        line([tvec0(1), tvec0(end)], [minspikeamp, minspikeamp], 'LineStyle', '--');    % line for minspikeamp
        text(.05, .15, ['Spike Frequency: ', num2str(spikefrequency), ' Hz'], 'Units', 'normalized');
        text(.05, .1, ['Spike Adaptation: ',  num2str(spikeadaptation), '%'], 'Units', 'normalized');
        text(.05, .05, ['Spike Threshold: ',  num2str(spikethreshold), 'mV'], 'Units', 'normalized');
    elseif ~isnan(ltspeaktime)            % LTS without bursts
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), vvec1(npi), 'bo', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr        % LTS by overrule
            plot(tvec0(npi), vvec1(npi), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                             % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                              % line for maxslope
        legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
            peakclass_label, peakfeature_label, 'Location', 'SouthOutside')
%        plot(tvec0(maxslopeind), vvec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tvec0(maxslopeind), vvec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakprom
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakwidth
    else                        % not LTS
        if peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tvec0(npi), vvec1(npi), 'kx', 'MarkerSize', 10);
        end
        if isempty(allspi)            % noise
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakclass_label, 'Location', 'SouthOutside')
        elseif isspontaneous            % spontaneous spikes
            plot(tvec0(allspi(1)), vvec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakclass_label, 'spontaneous spike', ...
                'Location', 'SouthOutside')
        else                    % spontaneous spike or noise
            plot(tvec0(allspi(1)), vvec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                peakclass_label, 'spontaneous spike or noise', ...
                'Location', 'SouthOutside')
        end
    end
    plot(tvec0(bl_start), vvec1(bl_start), 'g>');
    plot(tvec0(bl_end), vvec1(bl_end), 'y<');
    plot(tvec0(np_lbi_old), vvec3(np_lbi_old), 'k*');
    plot(tvec0(np_ubi_old), vvec3(np_ubi_old), 'k*');
    plot(tvec0(np_lbi), vvec3(np_lbi), 'r*');
    plot(tvec0(np_ubi), vvec3(np_ubi), 'r*');
    xlim([tvec0(1) tvec0(end)]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(strrep(filebase, '_', '\_'));
    figname = fullfile(outfolder, directories{1}, [filebase, '.png']);
    saveas(h, figname);

    % Plot scaled voltage trace for better comparison
    figure(h);
%    ylim([-120 20]);
    plot(tvec0, vvec1, 'g-', 'LineWidth', 1);
    xlimits = ltswin;
    xlim(xlimits);
    ylim([-100 -40]);            % Fix y-axis to determine whether the trace is good for fitting
    legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', peakclass_label);
    figname = fullfile(outfolder, directories{4}, [filebase, '_scaled.png']);
    saveas(h, figname);
    if peak2ndder > lts_thr_alt ...
        && peak2ndder <= lts_thr    % in "gray area"
        figname = fullfile(outfolder, directories{5}, [filebase, '_scaled.png']);
        saveas(h, figname);
    end
    if isempty(ptemp3) && ~isempty(ptemp2)    % LTS could be missed
        figname = fullfile(outfolder, directories{6}, [filebase, '_scaled.png']);
        saveas(h, figname);
        LTScouldbemissed = 1;
    end
    hold off;
%    close(h);

    % Plot LTS analysis
    % fprintf('Plotting LTS analysis ...\n');
    xlimits = ltswin;
    h = figure(5001);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'LTS analysis, moving-average-filtered trace');
    clf(h);
    subplot(3,1,1) % voltage trace
    plot(tvec0, vvec3, 'r-', 'LineWidth', 0.5); hold on
    if ~isnan(ltspeaktime) && ~isnan(bursttime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), vvec1(npi), 'go', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'gd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'ro', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr        % LTS by overrule
            plot(tvec0(npi), vvec1(npi), 'mo', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'md', 'MarkerSize', 8);
        end
    elseif ~isnan(ltspeaktime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), vvec1(npi), 'bo', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'bd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'ro', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr        % LTS by overrule
            plot(tvec0(npi), vvec1(npi), 'mo', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'md', 'MarkerSize', 8);
        end
    else
        if peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tvec0(npi), vvec1(npi), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tvec0(np_lbi), vvec3(np_lbi), 'r*');
    plot(tvec0(np_ubi), vvec3(np_ubi), 'r*');
    plot(tvec0(np_lbi_old), vvec3(np_lbi_old), 'k*');
    plot(tvec0(np_ubi_old), vvec3(np_ubi_old), 'k*');
    xlim(xlimits);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(['LTS analysis for ', strrep(filebase, '_', '\_'), ', moving-average-filtered trace']);
    subplot(3,1,2) % 1st derivative of voltage trace
    plot(tvec0(2:end), dvvec3, 'k-', 'LineWidth', 0.5); hold on
    plot(tvec0(2:end), dvvec3_sm, 'r-', 'LineWidth', 0.5);
    if ~isnan(ltspeaktime) && ~isnan(bursttime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), dvvec3(npi), 'go', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), dvvec3(maxslopeind), 'gd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), dvvec3(npi), 'ro', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), dvvec3(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr        % LTS by overrule
            plot(tvec0(npi), dvvec3(npi), 'mo', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), dvvec3(maxslopeind), 'md', 'MarkerSize', 8);
        end
    elseif ~isnan(ltspeaktime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), dvvec3(npi), 'bo', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), dvvec3(maxslopeind), 'bd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), dvvec3(npi), 'ro', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), dvvec3(maxslopeind), 'rd', 'MarkerSize', 8);
        elseif peak2ndder > lts_thr        % LTS by overrule
            plot(tvec0(npi), dvvec3(npi), 'mo', 'MarkerSize', 10);
            plot(tvec0(maxslopeind), dvvec3(maxslopeind), 'md', 'MarkerSize', 8);
        end
    else
        if peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), dvvec3(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tvec0(npi), dvvec3(npi), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tvec0(np_lbi), dvvec3(np_lbi), 'r*');
    plot(tvec0(np_ubi), dvvec3(np_ubi), 'r*');
    plot(tvec0(np_lbi_old), dvvec3(np_lbi_old), 'k*');
    plot(tvec0(np_ubi_old), dvvec3(np_ubi_old), 'k*');
    xlim(xlimits);
    xlabel('Time (ms)');
    ylabel('dV/dT');
    legend('unsmoothed', 'smoothed');
    subplot(3,1,3) % 2nd derivative of voltage trace
    plot(tvec0(3:end), ddvvec3, 'k-', 'LineWidth', 0.5); hold on
    line(xlimits, [lts_thr lts_thr], 'Color', 'r', ...
        'LineStyle', '--', 'LineWidth', 0.5);                       % mark LTS threshold
    line(xlimits, [lts_thr_alt lts_thr_alt], 'Color', 'b', ...
        'LineStyle', '--', 'LineWidth', 0.5);                       % mark alt LTS threshold
    xlim(xlimits);
    ax = gca;
    ylimits = get(ax, 'YLim');
    line([tvec0(ind3_begin) tvec0(ind3_begin)], ...
        ylimits, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);% mark ipeakt (start of LTS detection)
    if ~isnan(ltspeaktime) && ~isnan(bursttime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), ddvvec3(npi), 'go', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr        % in "gray area"
            plot(tvec0(npi), ddvvec3(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr         % LTS by overrule
            plot(tvec0(npi), ddvvec3(npi), 'mo', 'MarkerSize', 10);
        end
    elseif ~isnan(ltspeaktime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), ddvvec3(npi), 'bo', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr        % in "gray area"
            plot(tvec0(npi), ddvvec3(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr         % LTS by overrule
            plot(tvec0(npi), ddvvec3(npi), 'mo', 'MarkerSize', 10);
        end
    else
        if peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr        % in "gray area"
            plot(tvec0(npi), ddvvec3(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tvec0(npi), ddvvec3(npi), 'kx', 'MarkerSize', 10);
        end
    end
    plot(tvec0(np_lbi), ddvvec3(np_lbi), 'r*');
    plot(tvec0(np_ubi), ddvvec3(np_ubi), 'r*');
    plot(tvec0(np_lbi_old), ddvvec3(np_lbi_old), 'k*');
    plot(tvec0(np_ubi_old), ddvvec3(np_ubi_old), 'k*');
    xlabel('Time (ms)');
    ylabel('d2V/dT2');
    figname = fullfile(outfolder, directories{2}, [filebase, '_LTSanalysis', '.png']);
    saveas(h, figname);
    hold off;
    % close(h);

    % Plot burst analysis
    % fprintf('Plotting burst analysis ...\n');
    h = figure(5002);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Burst analysis, original trace');
    clf(h);
    plot(tvec0(np_lbi:np_ubi), vvec0(np_lbi:np_ubi), 'b-', 'LineWidth', 0.5); hold on
    plot(tvec2, vvec2, 'g.', 'LineWidth', 0.5);
    plot(tvec0(np_lbi:np_ubi), vvec3(np_lbi:np_ubi), 'r-', 'LineWidth', 0.5);
    if ~isnan(ltspeaktime) && ~isnan(bursttime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), vvec1(npi), 'go', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr        % LTS by overrule
            plot(tvec0(npi), vvec1(npi), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                             % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                              % line for maxslope
        plot(tvec0(bon_i), vvec0(bon_i), 'g>', 'MarkerSize', 10);       % triangle for burst onset time
        plot(tvec0(allspi), vvec0(allspi), 'gx', 'MarkerSize', 10);     % crosses for spikes
        legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
            peakclass_label, peakfeature_label, bursttime_label, spiketime_label, ...
            'Location', 'SouthOutside')
%        plot(tvec0(maxslopeind), vvec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tvec0(maxslopeind), vvec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakprom
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakwidth
        line([tvec0(1), tvec0(bon_i)], [spikethreshold, spikethreshold], ...
                'LineStyle', ':');                                      % line for spike threshold
        line([tvec0(1), tvec0(end)], [maxspikeamp, maxspikeamp], 'LineStyle', '--');    % line for maxspikeamp
        line([tvec0(1), tvec0(end)], [minspikeamp, minspikeamp], 'LineStyle', '--');    % line for minspikeamp
        text(.05, .7, ['Spike Frequency: ', num2str(spikefrequency), ' Hz'], 'Units', 'normalized');
        text(.05, .65, ['Spike Adaptation: ',  num2str(spikeadaptation), '%'], 'Units', 'normalized');
        text(.05, .6, ['Spike Threshold: ',  num2str(spikethreshold), 'mV'], 'Units', 'normalized');
    elseif ~isnan(ltspeaktime)
        if peak2ndder <= lts_thr_alt
            plot(tvec0(npi), vvec1(npi), 'bo', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'ro', 'MarkerSize', 10);
        elseif peak2ndder > lts_thr        % LTS by overrule
            plot(tvec0(npi), vvec1(npi), 'mo', 'MarkerSize', 10);
        end
%        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
%            'Color', 'b');                                % line for maxslope
        line([leftslopepoint_x, rightslopepoint_x], [leftslopepoint_y, rightslopepoint_y], ...
            'Color', 'm');                                % line for maxslope
        legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
            peakclass_label, peakfeature_label, 'Location', 'SouthOutside')
%        plot(tvec0(maxslopeind), vvec3(maxslopeind), 'm.');            % dot for maxslope
        plot(tvec0(maxslopeind), vvec1(maxslopeind), 'm.');             % dot for maxslope
        line([pkprom_x, pkprom_x], [pkprom_y1, pkprom_y2]);             % line for peakprom
        line([pkwidth_x1, pkwidth_x2], [halfprom_y, halfprom_y]);       % line for peakwidth
    else
        if peak2ndder > lts_thr_alt ...
            && peak2ndder <= lts_thr    % in "gray area"
            plot(tvec0(npi), vvec1(npi), 'rx', 'MarkerSize', 10);
        else
            plot(tvec0(npi), vvec1(npi), 'kx', 'MarkerSize', 10);
        end
        if isempty(allspi)
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakclass_label, 'Location', 'SouthOutside')
        elseif isspontaneous            % spontaneous spikes
            plot(tvec0(allspi(1)), vvec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakclass_label, 'spontaneous spike', ...
                'Location', 'SouthOutside')
        else                    % spontaneous spike or noise
            plot(tvec0(allspi(1)), vvec0(allspi(1)), 'rx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                peakclass_label, 'spontaneous spike or noise', ...
                'Location', 'SouthOutside')
        end
    end
    plot(tvec0(np_lbi_old), vvec3(np_lbi_old), 'k*');
    plot(tvec0(np_ubi_old), vvec3(np_ubi_old), 'k*');
    plot(tvec0(np_lbi), vvec3(np_lbi), 'r*');
    plot(tvec0(np_ubi), vvec3(np_ubi), 'r*');
%{
    ax = gca;
    if isspontaneous == 1
        text((19/20)*ax.XLim(1) + (1/20)*ax.XLim(2), ...
            (1/20)*ax.YLim(1) + (19/20)*ax.YLim(2), ...
            'Spontaneous spike!');
    end
%}
    xlim([tvec0(np_lbi), tvec0(np_ubi)]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(['Burst analysis for ', strrep(filebase, '_', '\_'), ', original trace']);
    figname = fullfile(outfolder, directories{3}, [filebase, '_burstanalysis', '.png']);
    saveas(h, figname);
    hold off;
    % close(h);
end

%{
%% OLD CODE

tvec2 = rsims*(1:round(tvec0(end)/rsims))';    % resampled time vector % NOT ROBUST

istart_ind = round(istart/sims);        % Index of current application % NOT ROBUST
ipeak_ind = round(ipeakt/sims);            % Index of current peak % NOT ROBUST

hr_ind = round(hrange(1)/sims):round(hrange(2)/sims);    % indices for calculating maxnoise % NOT ROBUST
ind3_begin = max(round(ltswin(1)/sims), ipeak_ind);    % first index for LTS search % NOT ROBUST
ind3_end = round(ltswin(2)/sims);            % last index for LTS search % NOT ROBUST

maxnoise = 2*std(vvec3(hr_ind));        % range of values of the median-filtered then 

plot(tvec0(maxslopeind), vvec1(maxslopeind), 'gd', 'MarkerSize', 8);    % max slope replaced by line

        [~, pkwidth_leftind] = min(abs(vvec3(1:npi)-(vvec3(npi)-peakprom/2)));

peakprom = vpeak_p3(psel3);            % this is the "minimum vertical distance that the signal must 
                        % descend on either side of the peak before either climbing back 
                        % to a level higher than the peak or reaching an endpoint" in mV
peakwidth = vpeak_w3(psel3) * sims;        % the width of the peak at half-prominence in ms
    ltspeakval = vvec1(npi);        % LTS amplitude, use median-filtered voltage trace

% Record indices relative to tvec0 or vvec0
npi = (v_begin - 1) + vpeak_i3(psel3);        % narrowest peak index
peaktime = (npi - istart_ind) * sims;        % narrowest peak time (delay)
% fprintf('Selected peak time == %g ms\n', peaktime);

    [maxslopeval, ind1] = max(dvvec3_sm(np_lbi:np_ubi));    % maximum slope in V/s, use moving-average-filtered voltage trace
    maxslopeind = np_lbi + ind1 - 1;            % index of maximum slope, relative to moving-average-filtered voltage trace

        leftslopepoint_y = maxslopeval * (tvec0(maxslopeind)-25 - tvec0(maxslopeind)) + vvec1(maxslopeind);
        rightslopepoint_y = maxslopeval * (tvec0(maxslopeind)+25 - tvec0(maxslopeind)) + vvec1(maxslopeind);

        line([tvec0(maxslopeind)-25, tvec0(maxslopeind)+25], [leftslopepoint_y, rightslopepoint_y]);

        leftslopepoint_x = (vvec3(np_ubi) - vvec1(maxslopeind)) / maxslopeval + tvec0(maxslopeind);
        rightslopepoint_x = (vvec1(npi) - vvec1(maxslopeind)) / maxslopeval + tvec0(maxslopeind);
        line([leftslopepoint_x, tvec0(maxslopeind)], [vvec3(np_ubi), vvec1(maxslopeind)]);
        line([tvec0(maxslopeind), rightslopepoint_x], [vvec1(maxslopeind), vvec1(npi)]);
% dv1 = diff(v1)/sims;
    % dv1_np_lbi = np_lbi - v_begin + 1;
    % dv1_np_ubi = np_ubi - v_begin + 1;
    % [maxslopeval, maxslopeind] = max(dv1(dv1_np_lbi:dv1_np_ubi));
    % maxslopeind = (maxslopeind + dv1_np_lbi - 1) + v_begin - 1;

    for x = 1+spacing_size:spacing_size:length(v1)
        if (v1(x)-v1(x-spacing_size))/(sims*spacing_size) > maxslopeval
            maxslopeval = (v1(x)-v1(x-spacing_size))/(sims*spacing_size);
            maxslopeind = x;
        end
    end
    maxslopeind = maxslopeind + v_begin - 1;

    % Find max slope based off of v1
    dv1 = zeros(length(v1), 1);
    maxslopeval = -realmax;
    maxslopeind = 1;
    spacing_size = 10;
    win_v1 = v1(np_lbi-v_begin+1:np_ubi-v_begin+1);
    for x = 1+spacing_size:length(win_v1)
        if (win_v1(x)-win_v1(x-spacing_size))/(sims*spacing_size) > maxslopeval
            maxslopeval = (win_v1(x)-win_v1(x-spacing_size))/(sims*spacing_size);
            maxslopeind = x;
        end
    end
    maxslopeind = maxslopeind + np_lbi;

            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'bd', 'MarkerSize', 8);
            plot(tvec0(maxslopeind), vvec1(maxslopeind), 'rd', 'MarkerSize', 8);

        leftslopepoint_x = (vvec3(np_ubi) - vvec1(maxslopeind)) / maxslopeval + tvec0(maxslopeind);
        rightslopepoint_x = (vvec1(npi) - vvec1(maxslopeind)) / maxslopeval + tvec0(maxslopeind);
        leftslopepoint_y = vvec3(np_ubi);
        rightslopepoint_y = vvec1(npi);
        line([leftslopepoint_x, rightslopepoint_x], [vvec3(np_ubi), vvec1(npi)], 'Color', 'm');    % line for maxslope

        leftslopepoint_x = tvec0(maxslopeind) - slope_spacing;
        rightslopepoint_x = tvec0(maxslopeind);
        leftslopepoint_y = vvec1(maxslopeind - spacing_size);
        rightslopepoint_y = vvec1(maxslopeind);

    % Moving-average-filter median-filtered traces for calculating maximum slope
    mafw3 = mafw3_dv/maxslopeval_appr;        % width in ms for the moving average filter for finding slopes
    ndp_mafw3 = round(mafw3/sims);
    vvec4 = smooth(vvec1, ndp_mafw3);
    v4 = vvec4(ind3);                % voltage vector of interest for calculating maxslope

if ~isempty(ptemp3) && isspontaneous ~= 1

% Burst onset index is the last point
%     lower than the base of the first spike
%     Find index in terms of v0
ind8 = find(v0(1:sp1sti(psel3)) < amp8, 1, 'last');

%}
