% m3ha_parse_dclamp_data.m
%% Parses all dynamic clamp data recorded by Christine
% Explanation:
%       Extract all dclamp data recorded by Christine between 
%           2010-09-17 to 2010-10-13 (11 folders' worth)
%       Reorganize data for each sweep and save sweep properties in a .mat file
%       Also save sweep properties in a .csv file 
%
% Requires:
%       cd/m3ha_append_lts_properties.m
%       cd/m3ha_find_ipsc_peak.m
%       cd/m3ha_find_ipsc_start_from_conductance.m
%       cd/m3ha_find_lts_many_sweeps.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_parse_mat.m
%       cd/m3ha_parse_sweep_settings.m
%       cd/find_passive_params.m
%       cd/m3ha_resave_sweeps.m
%       cd/m3ha_plot_correlations.m
%       cd/m3ha_plot_histograms_refine_threshold.m
%       cd/m3ha_compute_and_plot_statistics.m
%       cd/m3ha_estimate_passive_params.m

% File History:
% 2016-07-28 Modified from Christine's code
% 2016-09-09 Added flags, burst_thr, pk2sp_t, ptemp4, ptemp5
% 2016-09-09 Fixed fliplr in finding peak bounds to flipud
% 2016-09-09 Changed peak boundary detection into two steps; changed maxnoiseprom2 to 1 mV
% 2016-09-11 Changed maxnoiseprom2 to 2 mV and maxnoiseprom1 to 0.15 mV
% 2016-09-12 Changed maxnoiseprom1 to 0 mV
% 2016-09-13 Rendered IPSC offset obsolete by setting IPSC application to be 1000 ms for all sweeps.
% 2016-09-14 Added lts_thr_alt; changed minprom & maxnoiseprom2 to maxnoise (new statistic to save), 
%        which is equal to the baseline noise level (range of the mfmaf trace from 200~1000 ms)
% 2016-09-14 Changed trace for measuring holding potential from vvec0 to vvec1 
%        (due to spontaneous spikes messing up result, see E100810_0006_13 for example)
% 2016-09-14 Changed order of LTS determination: now it finds the first spike crossing 2nd derivative threshold
%        before checking whether it's a spontaneous spike
% 2016-09-14 Added peakclass as a statistic to save
% 2016-09-15 Changed alternative LTS threshold from -0.0065823 to -0.0086823
% 2016-09-15 Changed IPSC peak window from 1000~1300 ms to 1000~1350 ms
% 2016-09-15 Added peakprom as a statistic to save
% 2016-09-16 Changed IPSC peak window from 1000~1350 ms back to 1000~1300 ms
% 2016-09-16 Added vtraces_scaled
% 2016-09-22 figure(h) is needed before ylim in vtraces_scaled
% 2016-09-22 Fixed figure titles so that they include the trace name
% 2016-10-13 dclampdatalog_take4.csv now overwrites whenever saveswpinfoflag == 1
% 2016-10-13 Added maxslopetime & maxslopeval and plot them as diamonds
% 2016-10-13 Combined all the LTS detection to m3ha_find_lts.m under /home/Matlab/Adams_Functions
% 2016-10-13 m3ha_find_lts.m uses find(tvec0 >= istart, 1) instead of round(istart/sims), causing a slight difference in peaktimes
% 2016-10-13 Combined all the IPSC peak detection to m3ha_find_ipsc_peak.m under /home/Matlab/Adams_Functions
% 2016-10-14 Added peakwidth as a statistic to save
% 2016-10-15 Added functionsdirectory, homedirectory and addpath
% 2016-10-18 Added close all inside the parfor loop
% 2016-10-18 Made functionsdirectory, homedirectory dependent on existence
% 2016-10-27 BT - Added 'maxspikeamp', 'minspikeamp', 'spikefrequency', 'spikeadaptation'
% 2016-10-27 Combined all the Rinput detection to find_passive_params.m under /home/Matlab/Adams_Functions
% 2016-10-31 Renamed narrowpeaktime & narrowpeak2ndder as peaktime & peak2ndder
% 2016-10-31 Added m3ha_append_lts_properties.m, which generates vectors of peak features restricted to those with LTS
% 2016-11-01 Added compute flags
% 2016-11-01 Added plotpassive2flag (for m3ha_estimate_passive_params.m)
% 2016-11-07 Added cpa_ap, g_sc, i_sc
% 2016-11-07 Reorganized code to make more efficient; created m3ha_parse_sweep_settings.m & m3ha_resave_sweeps.m
% 2016-11-30 maxsets was changed from 347 to 346; maxswps was changed from 7455 to 7430
% 2016-12-13 Reversed sign of LJP
% 2017-01-16 Changed current pulse response to last just 150 ms (cprwin is changed from [95, 500] to [95, 260])
% 2017-02-16 Changed peakclass_label to nine items
% 2017-03-22 BT - Added LTScouldbemissed
% 2017-12-21 Change tabs to spaces
% 2017-12-21 Added spikethreshold, firstspiketime and lastspiketime
% 2018-08-06 Made /tmp/data/m3ha/ the first priority home directory
% 2018-08-06 Updated usage of m3ha_parse_sweep_settings.m
% 2018-10-04 Now uses m3ha_parse_mat.m

%% Flags
debugflag = 0;
resavedataflag = 0;
computepassiveflag = 1;
plotpassiveflag = 0;
computeIPSCoffsetoldflag = 1;
plotIPSCoffsetoldflag = 0;
computeIPSCpeakflag = 1;
plotIPSCpeakflag = 0;
computeLTSflag = 1;
plotLTSflag = 0; %1;
saveswpinfoflag = 1;
plotpassive2flag = 0; %1;
plothistogramsflag = 0; %1;
plotcorrelationsflag = 0; %1;
plotbargraphsflag = 0; %1;
preallocateflag = 1;    % Only use this after total number of sweeps and total number of cells are certain
                        % Can't use this under debug mode

%% Fixed parameters used in the experiments
ljp = 10;               % Liquid junction potential used in Christine's experiments (mV)

cpwin = [95, 115];      % Window in which the current pulse would lie (ms) 
                        % (Supposed to be 100-110 ms but there will be offset)
cprwin = [95, 260];     % Window in which the current pulse response would lie (ms)
cpmid = 105;            % Approximate midpoint of the current pulse (ms)
ipsctwin = [1000, 1100];% Window in which IPSC is first applied (ms) (obsolete)
IPSC_start_time = 1000; % Time at which IPSC is first applied (ms)
ipscpwin = [1000, 1300];% Window in which IPSC reaches peak amplitude (ms) 
                        %     Based on observation, IPSCs are not influenced by LTSs before 1300 ms
hrange = [200, 1000];   % Range in which holding current should maintain at a certain holding potential
ltswin = [1000, 7960];  % Window in which the low threshold spike would lie (ms) 
                % Previously [1000 4500] by Christine;
                % 8000 ms is the approximate time that IPSC is terminated; subtract out 40 ms for discontinuities
                % 1000 ms will be replaced by ipscpt, the time when IPSC reaches peak amplitude
if preallocateflag
    maxcells = 49;      % Total number of neurons recorded
    maxsets = 346;      % Total number of cell-pharm-vhold conditions recorded
    maxswps = 7430;     % Total number of sweeps recorded
end

%% Set sweep properties to log
logheader = {'Data filename', 'Cell ID #', ...
        'Pharm condition', 'Vhold (mV)', 'GABAB IPSC G incr (%)', 'Within condition sweep #', ...
        'GABAB IPSC G amp (nS)', 'GABAB IPSC G Trise (ms)', 'GABAB IPSC G TfallFast (ms)', ...
        'GABAB IPSC G TfallSlow (ms)', 'GABAB IPSC G w', ...
        'Current pulse amplitude (pA)', 'Rinput (MOhm)', ...
        'IPSC offset (obsolete) (ms)', 'IPSC peak time (ms)', 'IPSC peak amplitude (pA)', ...
        'Actual Ihold (pA)', 'Actual Vhold (mV)', 'Maximum noise (mV)', ...
        'Peak time (ms)', 'Peak 2nd derivative (V^2/s^2)', ...
        'Peak prominence (mV)', 'Peak width at half-prom (ms)', 'Peak class #', 'Spikes per peak', ...
        'LTS onset time (ms)', 'LTS amplitude (mV)', ...
        'LTS maximum slope time (ms)', 'LTS maximum slope amplitude (V/s)', ...
        'Burst onset time (ms)', 'Spikes per burst', ...
        'Spike threshold (mV)', 'First spike time (ms)', 'Last spike time (ms)', ...
        'Maximum spike amplitude (mV)', 'Minimum spike amplitude (mV)', ...
        'Spike frequency (Hz)', 'Spike adaptation (%)'};
logvariables = {'fnrow', 'cellidrow', ...
        'prow', 'vrow', 'grow', 'swpnrow', ...
        'gabab_amp', 'gabab_Trise', 'gabab_TfallFast', ...
        'gabab_TfallSlow', 'gabab_w', ...
        'currpulse', 'Rin', 'ioffset_old', 'imint', 'imin', ...
        'actIhold', 'actVhold', 'maxnoise', 'peaktime', 'peak2ndder', ...
        'peakprom', 'peakwidth', 'peakclass', 'spikesperpeak', ...
        'ltspeaktime', 'ltspeakval', 'maxslopetime', 'maxslopeval', ...
        'bursttime', 'spikesperburst', ...
        'spikethreshold', 'firstspiketime', 'lastspiketime', ...
        'maxspikeamp', 'minspikeamp', ...
        'spikefrequency', 'spikeadaptation'};

%% Peak Classification (corresponds to peakclass == 1 ~ 9):
%    Note: this variable is not used here but a modified version is applied in m3ha_find_lts.m
peakclassLabels = {'Not LTS by prominence', ...
            'Not LTS by narrowness', ...
            'Not LTS by shape', ...
            'Not LTS by overrule', ...
            'LTS with no burst by overrule', ...
            'LTS with no burst; contentious', ...
            'LTS with burst; contentious', ...
            'LTS with no burst; definite', ...
            'LTS with burst; definite'};


%% Files to test for debug mode
%{
% Previous problematic files
debug_files = {'D091710_0007', 'A092110_0005', 'A092110_0012', 'A092110_0013', ...
        'A092810_0001', 'A092810_0002', 'A092810_0004', ...
        'A092810_0005', 'A092810_0006', 'A092810_0007', ...
        'A092910_0001', 'A092910_0002', 'B092710_0007', ...
        'D091710_0013', 'G091810_0003', 'A092110_0012', ...
        'B092710_0007', 'A092110_0005', 'E100810_0003', ...
        'D092910_0007'};
%}
% 2016-09-09 Fixed Not_the_first_LTS by 
%        finding the first LTS peak that crosses 2nd derivative threshold
%         instead of the narrowest peak
% debug_files = {'E101210_0005', 'G101210_0006'};

% 2016-09-09 Fixed Missed_LTS by making exceptions in finding ptemp3
% debug_files = {'A092910_0001', 'F101210_0000'};

% 2016-09-11 Fixed Spikes_per_burst_incorrect by changing maxnoiseprom2 to 2 mV
% debug_files = {'G101210_0001', 'G101210_0003', 'G101210_0006'};

% 2016-09-12 Fixed Spontaneous_spikes_too_close_together by fixing peak boundary detection (changing fliplr to flipud)
%         & changing maxnoiseprom1 to 0 mV
%         Note: When maxnoiseprom1 was allowed to be 0.15, only one file in this condition was problematic
% debug_files = {'A100110_0008', 'A100110_0011', 'A100810_0000', ...
%        'C092710_0001', 'D092710_0003', 'D092710_0005', ...
%        'D092710_0006', 'D101310_0008', 'E100810_0000', ...
%        'E100810_0001', 'E100810_0002', 'E100810_0003', ...
%        'E100810_0004', 'E100810_0006', 'F092910_0001', ...
%        'J101210_0000', 'J101210_0003'};

% LTS peak before first action potential
% debug_files = {'A092910_0001', 'F101210_0000', 'G101210_0001', 'G101210_0001'};

% Small peak before LTS peak
% debug_files = {'B101210_0009', 'F092710_0006'};

% debug_files = {'E092910_0000', 'J101210_0000'};

% Wrong holding potential
% debug_files = {'E100810_0006'};

% debug_files = {'A092910_0001', 'F101210_0000', 'G101210_0001', 'G101210_0001', ...
%        'B101210_0009', 'F092710_0006', 'A092810_0004', ...
%        'A092810_0005', 'A092810_0006', 'A092810_0007', ...
%        'A092910_0001', 'A092910_0002'};

% debug_files = {'E091710_0000'};

debug_files = {'A092810_0006', 'E091710_0000'};

%{
debug_files = {'F092710_0006'};                                 % Missed_LTS_by_order
debug_files = {'G101210_0001', 'F101210_0000', 'A092910_0001'}; % Missed_LTS_by_shape
debug_files = {'E101210_0002'};                                 % Spikes_per_burst_incorrect
debug_files = {'D092710_0007', 'F092710_0006', 'F092710_0007'}; % Noise_in_trace
debug_files = {'A092910_0008', 'D092710_0005', 'E100110_0004'}; % Spontaneous_LTSs_or_bursts
debug_files = {'A100810_0000', 'E092710_0000', 'G101310_0000'}; % Wide_LTS_could_be_noise
debug_files = {'A092110_0005', 'A101210_0000', 'J101210_0004'}; % Looks_like_LTS_not_by_narrowness
debug_files = {'D101310_0008', 'F092710_0010', 'G101310_0000'}; % Looks_like_LTS_not_by_prominence
debug_files = {'C092910_0001', 'E100810_0000'};                 % Looks_like_missed_LTS
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make sure flags are consistent
if debugflag && preallocateflag
    preallocateflag = 0;
    fprintf('WARNING: preallocateflag cannot be used in debug mode\n')
end

%% Locate home directory
homeDirectory = m3ha_locate_homedir;

%% Set folders for reading and saving files
infolder = fullfile(homedirectory, '/data_dclamp/mark_dclamp/');
newinfolder = fullfile(homedirectory, '/data_dclamp/take4/matfiles/');
if debugflag
    outfolder = fullfile(homedirectory, '/data_dclamp/take4/debug/');
else
    outfolder = fullfile(homedirectory, '/data_dclamp/take4/');
end
check_dir(outfolder)

%% Count total number of usable cells, sets and sweeps, generate filenames and record sweep properties
if preallocateflag == 1
    [numcells, numsets, numswps, ...
        cellnames, abffullfn, nswps_cpv, nswps_used, ...
        fnrow, cellidrow, prow, vrow, grow, swpnrow, ...
        gabab_amp, gabab_Trise, gabab_TfallFast, gabab_TfallSlow, gabab_w, actIhold] = ...
        m3ha_parse_sweep_settings (infolder, ljp, debugflag, debug_files, preallocateflag, maxcells, maxsets, maxswps);
else
    [numcells, numsets, numswps, ...
        cellnames, abffullfn, nswps_cpv, nswps_used, ...
        fnrow, cellidrow, prow, vrow, grow, swpnrow, ...
        gabab_amp, gabab_Trise, gabab_TfallFast, gabab_TfallSlow, gabab_w, actIhold] = ...
        m3ha_parse_sweep_settings (infolder, ljp, debugflag, debug_files, preallocateflag);
end

%% Extract .abf data and resave as .mat file
if resavedataflag == 1
    % Create output directory if not already exist
    check_subdir(outfolder, 'matfiles');

    % Initialize vectors
    cpa_ap = cell(1, numsets);
    g_sc = cell(1, numsets);
    i_sc = cell(1, numsets);

    % FOR each set (each cell-pharm-Vhold condition)
    for sn = 1:numsets

        % Start timer
        tic;

        % Extract name of set and number of sweeps used
        [~, datfn, ~] = fileparts(abffullfn{sn});
        nswps = nswps_cpv(sn);      % total number of sweeps in this set
        n_used = nswps_used(sn);    % total number of sweeps used before this set

        % Find GABAB IPSC parameters for this set
        gabab_amp_now = gabab_amp(n_used+1:n_used+nswps);
        gabab_Trise_now = gabab_Trise(n_used+1:n_used+nswps);
        gabab_TfallFast_now = gabab_TfallFast(n_used+1:n_used+nswps);
        gabab_TfallSlow_now = gabab_TfallSlow(n_used+1:n_used+nswps);
        gabab_w_now = gabab_w(n_used+1:n_used+nswps);

        % Extract .abf data and resave as .mat file
        fprintf(['RESAVING .abf into .mat files for the set ', datfn, ' ...\n']);
        [cpa_ap{sn}, g_sc{sn}, i_sc{sn}] = m3ha_resave_sweeps(abffullfn{sn}, nswps, ljp, ...
            IPSC_start_time, gabab_amp_now, gabab_Trise_now, ...
            gabab_TfallFast_now, gabab_TfallSlow_now, gabab_w_now, cpmid, outfolder);

        % End timer
        toc;
        fprintf('\n');
    end
end


%% Analyze data and record statistics
% Create output directories if not already exist
if plotpassiveflag == 1
    check_subdir(outfolder, 'passive');

end
if plotIPSCoffsetoldflag == 1
    check_subdir(outfolder, 'IPSCoffset_old');
end
if plotIPSCpeakflag == 1
    check_subdir(outfolder, 'IPSCpeak');
end
if plotLTSflag == 1
    check_subdir(outfolder, ...
                {'vtraces', 'LTSanalysis', 'burstanalysis', ...
                    'vtraces_scaled', 'gray_area', 'LTScouldbemissed'});
end

% Initialize variables different for each set
swp_ind = cell(1, numsets);     % sweep indices for each set
cpa_mean = cell(1, numsets);    % mean current pulse amplitude (nA)
dvss = cell(1, numsets);        % steady state voltage difference after current step (mV)
IPSC_offset_old = cell(1, numsets); % IPSC offset using conductance trace (obsolete)

cpa = cell(1, numsets);         % current pulse amplitude (pA)
Rinput = cell(1, numsets);      % input resistance (MOhm)
IPSC_offset_old2 = cell(1, numsets);% IPSC offset using current trace (obsolete)
ipeakdelay = cell(1, numsets);  % IPSC peak delay (ms)
ipeak_amp = cell(1, numsets);   % IPSC peak amplitude (pA)
bl_mean = cell(1, numsets);     % Actual Vhold (mV)
mnoise = cell(1, numsets);      % Maximum noise (mV)
npt = cell(1, numsets);         % Narrowest peak time (ms)
np2der = cell(1, numsets);      % Narrowest peak 2nd derivative (V^2/s^2)
pk_prom = cell(1, numsets);     % Peak prominence (mV)
pk_width = cell(1, numsets);    % Peak width (ms)
pk_class = cell(1, numsets);    % Peak class #
spp = cell(1, numsets);         % Spikes per peak
ltst = cell(1, numsets);        % LTS peak time (ms)
ltsv = cell(1, numsets);        % LTS peak value (mV)
mxslt = cell(1, numsets);       % LTS maximum slope time (ms)
mxslv = cell(1, numsets);       % LTS maximum slope value (mV)
btime = cell(1, numsets);       % Burst onset time (ms)
spb = cell(1, numsets);         % Spikes per burst
maxspi = cell(1, numsets);      % Maximum spike amplitude (mV)
minspi = cell(1, numsets);      % Minimum spike amplitude (mV)
spif = cell(1, numsets);        % Spikes frequency (Hz)
spia = cell(1, numsets);        % Spike adaptation (%)

% FOR each set (each pharm-Vhold pair)
parfor sn = 1:numsets
%for sn = 1:numsets
    % Start timer
    tic;

    % Extract name of set and number of sweeps
    [~, datfn, ~] = fileparts(abffullfn{sn});
    nswps = nswps_cpv(sn);      % total number of sweeps in this set
    n_used = nswps_used(sn);    % total number of sweeps used before this set
    ind_first = n_used + 1;     % index of first sweep
    ind_last = n_used + nswps;  % index of last sweep
    swp_ind{sn} = ind_first:ind_last;   % sweep indices for this set

    % Load relevant vectors from mat files
    fprintf(['loading .mat files for the set ', datfn, ' ...\n']);

    % Get all the matfiles for this set
    matfiles = arrayfun(@(x) fullfile(newinfolder, ...
                                [datfn, '_', num2str(x), '.mat']), ...
                        1:nswps, 'UniformOutput', false);

    % Parse and load matfiles
    [parsedParams, parsedData] = m3ha_parse_mat(matfiles);
    if isempty(parsedParams)
        error('Something wrong with the set %s', datfn);
    end

    % Extract parameters
    ndps = parsedParams.ndps;
    sims = parsedParams.sims;
    ndps2 = parsedParams.ndps2;
    sims2 = parsedParams.sims2;

    % Extract vectors
    tvec0 = parsedData.tvec0;
    tvec2 = parsedData.tvec2;
    gvec0s = parsedData.gvec0s;
    ivec0s = parsedData.ivec0s;
    vvec0s = parsedData.vvec0s;
    gvec1s = parsedData.gvec1s;
    ivec1s = parsedData.ivec1s;
    vvec1s = parsedData.vvec1s;
    gvec2s = parsedData.gvec2s;
    ivec2s = parsedData.ivec2s;
    vvec2s = parsedData.vvec2s;
    vvec3s = parsedData.vvec3s;

    if computepassiveflag
        % Analyze passive parameters such as input resistance (MOhm)
        fprintf('ANALYZING passive parameters for %s ...\n', datfn);
        [params, cpa{sn}] = ...
            find_passive_params (tvec0, ivec0, vvec0, ...
                'PulseWindow', cpwin, 'PulseResponseWindow', cprwin, ...
                'PlotFlag', plotpassiveflag, 'OutFolder', outfolder, ...
                'FileBase', datfn, 'Ivec1s', ivec1);

        cpa{sn} = params.pulseAmplitude;
        cpa_mean{sn} = params.cpa_mean * ones(1, nswps);% mean current pulse amplitude (nA)
        dvss{sn} = params.dvss * ones(1, nswps);        % steady state voltage difference after current step (mV)
        Rinput{sn} = params.Rinput * ones(1, nswps);    % input resistance (MOhm)
    end
    if computeIPSCoffsetoldflag
        % Analyze time of IPSC start (obsolete)
        fprintf('FINDING (and plotting) IPSC offsets (obsolete) ...\n');
        [IPSC_offset1, IPSC_offset2, ~, ~] = ...
            m3ha_find_ipsc_start_from_conductance (tvec0, gvec1, ivec0, ipsctwin, plotIPSCoffsetoldflag, vvec0, outfolder, datfn);

        IPSC_offset_old{sn} = IPSC_offset1 * ones(1, nswps);     
        IPSC_offset_old2{sn} = IPSC_offset2 * ones(1, nswps);    
    end
    if computeIPSCpeakflag
        % Find and plot IPSC peaks
        fprintf('FINDING (and plotting) IPSC peaks ...\n');
        [ipeak_time, ~, ipeak_amp{sn}, ipeakdelay{sn}] = ...
            m3ha_find_ipsc_peak(tvec0, ivec0, IPSC_start_time, ...
                ipscpwin, plotIPSCpeakflag, outfolder, datfn);

    end
    if computeLTSflag
        % Find and plot LTSs
        fprintf('FINDING (and plotting) LTSs ...\n');
        [bl_mean{sn}, mnoise{sn}, npt{sn}, np2der{sn}, pk_prom{sn}, pk_width{sn}, pk_class{sn}, ...
        spp{sn}, ltsv{sn}, ltst{sn}, mxslv{sn}, mxslt{sn}, btime{sn}, spb{sn}, ...
        spthr{sn}, fspt{sn}, lspt{sn}, maxspi{sn}, minspi{sn}, ...
        spif{sn}, spia{sn}] = ...
            m3ha_find_lts_many_sweeps (tvec0, vvec0, IPSC_start_time, ...
                ipeak_time, hrange, ltswin, ...
                plotLTSflag, outfolder, datfn, ...
                tvec2, vvec1, vvec2, vvec3);
    end

    % Stop timer
    toc;
    fprintf('\n');

    close all;
end

% Record variables different for each sweep
if resavedataflag
    currpulse_appr = cell2mat(cpa_ap);
    condscale = cell2mat(g_sc);
    currscale = cell2mat(i_sc);

    % Save variables
    filename = 'dclampdatalog_take4_resave.mat';
    resavedatafn = fullfile(outfolder, filename);
    save(resavedatafn, 'currpulse_appr', 'condscale', 'currscale', '-v7.3');
end
if computepassiveflag
    currpulse = cell2mat(cpa);          % Current pulse amplitude (pA)
    Rin = cell2mat(Rinput);             % Rinput (MOhm)
end
if computeIPSCoffsetoldflag
    ioffset_old = cell2mat(IPSC_offset_old2);    % IPSC offset (ms)
end
if computeIPSCpeakflag
    imint = cell2mat(ipeakdelay);       % IPSC peak delay (ms)
    imin = cell2mat(ipeak_amp);         % IPSC peak amplitude (pA)
end
if computeLTSflag
    actVhold = cell2mat(bl_mean);       % Actual Vhold (mV)
    maxnoise = cell2mat(mnoise);        % Maximum noise (mV)
    peaktime = cell2mat(npt);           % Narrowest peak time (ms)
    peak2ndder = cell2mat(np2der);      % Narrowest peak 2nd derivative (V^2/s^2)
    peakprom = cell2mat(pk_prom);       % Peak prominence (mV)
    peakwidth = cell2mat(pk_width);     % Peak width (ms)
    peakclass = cell2mat(pk_class);     % Peak class #
    spikesperpeak = cell2mat(spp);      % Spikes per peak
    ltspeaktime = cell2mat(ltst);       % LTS peak time (ms)
    ltspeakval = cell2mat(ltsv);        % LTS peak value (mV)
    maxslopetime = cell2mat(mxslt);     % LTS maximum slope time (ms)
    maxslopeval = cell2mat(mxslv);      % LTS maximum slope value (mV)
    bursttime = cell2mat(btime);        % Burst onset time (ms)
    spikesperburst = cell2mat(spb);     % Spikes per burst
    spikethreshold = cell2mat(spthr);   % Spike threshold (mV)
    firstspiketime = cell2mat(fspt);    % First spike time (ms)
    lastspiketime = cell2mat(lspt);     % Last spike time (ms)
    maxspikeamp = cell2mat(maxspi);     % Maximum spike amplitude (mV)
    minspikeamp = cell2mat(minspi);     % Minimum spike amplitude (mV)
    spikefrequency = cell2mat(spif);    % Spikes frequency (Hz)
    spikeadaptation = cell2mat(spia);   % Spike adaptation (%)
end

%% Save variables that store sweep properties
% Must be consistent with logheader & logvariables
filename = 'dclampdatalog_take4.mat';
swpdatafn = fullfile(outfolder, filename);
if saveswpinfoflag == 1
    command = ['save(swpdatafn, ', ...
        '''logheader'', ''logvariables'', ', ...
        '''cellnames'', ''abffullfn'', ''nswps_cpv'', ''nswps_used'', '];
    for f = 1:numel(logvariables)
        command = [command, sprintf('''%s'', ', logvariables{f})];
    end
    command = [command, '''-v7.3'');'];
    eval(command);

    % Print sweep properties into a comma-separated-value datalog file
    % Must be consistent with logheader & logvariables
    logpath = fullfile(outfolder, 'dclampdatalog_take4.csv');
    fid = fopen(logpath, 'w');
    for k = 1:numel(logheader)
        if k < numel(logheader)
            fprintf(fid, '%s, ', logheader{k});
        else
            fprintf(fid, '%s\n', logheader{k});
        end
    end
    for k = 1:numel(logvariables)
        if k < numel(logvariables)
            fprintf(fid, '%s, ', logvariables{k});
        else
            fprintf(fid, '%s\n', logvariables{k});
        end
    end
    for k = 1:numswps
        fprintf(fid, ['%s, %d, ' ...
        '%d, %d, %d, %d, ' ...
        '%g, %g, %g, ', ...
        '%g, %g, ' ...
        '%g, %g, %g, %g, %g, ' ...
        '%g, %g, %g, %g, %g, ' ...
        '%g, %g, %d, %d, ' ...
        '%g, %g, %g, %g, ' ...
        '%g, %d, ', ...
        '%g, %g, %g, ', ...
        '%g, %g, ', ...
        '%g, %g\n'], ...
        fnrow{k}, cellidrow(k), ...
        prow(k), vrow(k), grow(k), swpnrow(k), ...
        gabab_amp(k), gabab_Trise(k), gabab_TfallFast(k), ...
        gabab_TfallSlow(k), gabab_w(k), ...
        currpulse(k), Rin(k), ioffset_old(k), imint(k), imin(k), ...
        actIhold(k), actVhold(k), maxnoise(k), peaktime(k), peak2ndder(k), ...
        peakprom(k), peakwidth(k), peakclass(k), spikesperpeak(k), ...
        ltspeaktime(k), ltspeakval(k), maxslopetime(k), maxslopeval(k), ...
        bursttime(k), spikesperburst(k), ...
        spikethreshold(k), firstspiketime(k), lastspiketime(k), ...
        maxspikeamp(k), minspikeamp(k), ...
        spikefrequency(k), spikeadaptation(k));
    end
    fclose(fid);
end

%% Generate vectors of peak features restricted to those with LTS
m3ha_append_lts_properties(swpdatafn);

%% Refine spike threshold
noburstind = find(isnan(bursttime));        % indices of sweeps with no bursts
max_lts_amp = max(ltspeakval(noburstind));  % maximum LTS amplitude without bursts 
fprintf('The highest LTS peak without bursts has amplitude == %g mV\n', max_lts_amp);
fprintf('The files with such large LTS peak amplitude are:\n');
fnrow(noburstind(find(ltspeakval(noburstind) >= -48.1)))

%% Extract passive parameters for each cell, Vhold set
if plotpassive2flag == 1
    m3ha_estimate_passive_params(0, outfolder, outfolder);
    m3ha_estimate_passive_params(1, outfolder, outfolder);
    if ~debugflag
        m3ha_estimate_passive_params(2, outfolder, outfolder);
    end
end

%% Plot histograms and refine threshold to use to define an LTS
if plothistogramsflag == 1
    fprintf('Plotting histograms and refining threshold ...\n');
    m3ha_plot_histograms_refine_threshold(0, outfolder, outfolder);
    m3ha_plot_histograms_refine_threshold(1, outfolder, outfolder);
    if ~debugflag
        m3ha_plot_histograms_refine_threshold(2, outfolder, outfolder);
    end
end

%% Plot correlation diagrams
if plotcorrelationsflag == 1
    fprintf('Plotting correlation diagrams ...\n');
    m3ha_plot_correlations(0, outfolder, outfolder);
    m3ha_plot_correlations(1, outfolder, outfolder);
    if ~debugflag
        m3ha_plot_correlations(2, outfolder, outfolder);
    end
end

%% Plot bar graphs
if plotbargraphsflag == 1
    fprintf('Plotting bar graphs ...\n');
    m3ha_compute_and_plot_statistics('DataMode', 0, 'Directory', outfolder);
    m3ha_compute_and_plot_statistics('DataMode', 1, 'Directory', outfolder);
    if ~debugflag
        m3ha_compute_and_plot_statistics('DataMode', 2, 'Directory', outfolder);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
%% OLD CODE

        % The function m3ha_find_ipsc_peak.m is under /home/Matlab/Adams_Functions
        % The function m3ha_find_lts.m is under /home/Matlab/Adams_Functions

    % Record analyzed sweep properties
    fprintf('Recording analyzed sweep properties ...\n');

function dclampDataExtractor (varargin)

cpwin = [95, 115];        % Window in which the current pulse would lie (ms) (Supposed to be 100-110 ms but there will be offset)


        % Record analyzed sweep properties
        for swp = 1:nswps        % FOR each sweep
            swpc = n_used + swp;            % overall sweep count;
        end

plotRinputflag = 1;
                    first_dip_pt = cell(1, nswps);
                    cpa = zeros(1, nswps);
                    dv0 = zeros(1, nswps);
                    cfit = cell(1, nswps);
                    xvecs = cell(1, nswps);
                    yvecs = cell(1, nswps);
                    dv1 = zeros(1, nswps);
                    ind = round(cpwin(1)/sims):round(cpwin(2)/sims);% indices of interest
                    ind1 = round(cpwin(1)/sims):round(cpmid/sims);    % indices for finding current pulse start
                    ind2 = round(cpmid/sims):round(cpwin(2)/sims);    % indices for finding current pulse end
                    mvind = 1:round(mvw/sims);            % indices for taking the mean of voltages
                    base_ind = zeros(nswps, length(mvind));
                    last_ind = zeros(nswps, length(mvind));
                    ivec1_part1_begin = ind1(1);
                    ivec1_part2_begin = ind2(1);
                    ivec1_part1 = ivec1(ind1, :);            % Use median-filtered current trace
                    ivec1_part2 = ivec1(ind2, :);
                    ivec1_part3 = ivec1(round(cpwin(1)/sims):round((cpwin(1) + mvw)/sims), :);
                                            % indices for measuring cp baseline
                    ivec1_part4 = ivec1(round(cpmid/sims):round((cpmid + mvw)/sims), :);
                                            % indices for measuring cp peak
                    ft = fittype('a*exp(-x/b)+c*exp(-x/d)+e');    % double exponential
                    coeff = coeffnames(ft);                % {'a'; 'b'; 'c'; 'd'; 'e'}
                    parfor swp = 1:nswps                % FOR each sweep
                        cpa(swp) = mean(ivec1_part4(:,swp)) - mean(ivec1_part3(:,swp));    % current pulse amplitude (pA)
                        first_dip_pt{swp} = find(ivec1_part1(:, swp) > cpa(swp)/4, 1, 'last');
                        if isempty(first_dip_pt{swp})
                            xvecs{swp} = [];
                            yvecs{swp} = [];
                        else
                            % Find the voltage trace corresponding to the current pulse
                            cpstart = (ivec1_part1_begin - 1) + first_dip_pt{swp};    % index of current pulse start
                            before_rise_pt = find(ivec1_part2(:, swp) < cpa(swp) * 3/4, 1, 'last');
                            cpend = (ivec1_part2_begin - 1) + before_rise_pt;    % index of current pulse end
                            base_ind(swp, :) = cpstart - round(0.5/sims) - fliplr(mvind);    % base indices
                            last_ind(swp, :) = cpend - fliplr(mvind);    % last indices of the current pulse
                            basev = mean(vvec0(base_ind(swp, :), swp));
                            lastv = mean(vvec0(last_ind(swp, :), swp));
                            dv0(swp) = basev - lastv; % change in membrane potential on voltage trace (mV)
                            if dv0(swp) <= 0
                                xvecs{swp} = [];
                                yvecs{swp} = [];
                            else
                                xvecs{swp} = tvec0(cpstart:cpend) - tvec0(cpstart);
                                yvecs{swp} = vvec0(cpstart:cpend, swp) - basev;
                            end
                        end
                    end
                    dv0_mean = mean(dv0(dv0 > 0));
                    cpa_mean(sn) = mean(cpa(cpa < 0));
                    xvec_all = [];    % put data from all sweeps together
                    yvec_all = [];
                    for swp = 1:nswps
                        xvec_all = cat(1, xvec_all, xvecs{swp});
                        yvec_all = cat(1, yvec_all, yvecs{swp});
                    end
                    % Fit voltage trace with two exponentials
                    if isempty(xvec_all)
                        dv1(sn) = NaN;
                        Rinput(sn) = NaN;        % To be detected in .csv filename
                    else
                        cfit = fit(xvec_all, yvec_all, ft, ...
                            'StartPoint', [dv0_mean, 2, dv0_mean, 2, -2 * dv0_mean], ...
                            'Lower', [0, 0, 0, 0, -20 * dv0_mean], ...
                            'Upper', [10 * dv0_mean, 100, 10 * dv0_mean, 100, 0]); 
                            % typical tau ~ 2 ms
                        dv1(sn) = -(cfit.a + cfit.c);
                        Rinput(sn) = (dv1(sn)*10^-3)/(cpa_mean(sn)*10^-12)/10^6;    
                            % input resistance in MOhm
                    end

                    if plotRinputflag == 1
                        h = figure('Visible', 'off');
                        set(h, 'Name', 'Input resistance analysis');
                        clf(h);
                        subplot(3,1,1);
                        for swp = 1:nswps            % Plot each current trace
                            plot(tvec0(ind), ivec0(ind, swp), 'b'); hold on; 
                            plot(tvec0(ind), ivec1(ind, swp), 'g'); 
                            if ~isempty(first_dip_pt{swp})
                                plot(tvec0(base_ind(swp, 1)), ivec0(base_ind(swp, 1), swp), '>');
                                plot(tvec0(base_ind(swp, end)), ivec0(base_ind(swp, end), swp), '<');
                                plot(tvec0(last_ind(swp, 1)), ivec0(last_ind(swp, 1), swp), '>');
                                plot(tvec0(last_ind(swp, end)), ivec0(last_ind(swp, end), swp), '<');
                            end
                        end
                        title(['Input resistance analysis for ', datfn])
                        xlabel('Time (ms)')
                        ylabel('Current (pA)')
                        xlim(cpwin);
                        subplot(3,1,2);
                        for swp = 1:nswps            % Plot each voltage trace
                            plot(tvec0(ind), vvec0(ind, swp)); hold on; 
                            if ~isempty(first_dip_pt{swp})
                                plot(tvec0(base_ind(swp, 1)), vvec0(base_ind(swp, 1), swp), '>');
                                plot(tvec0(base_ind(swp, end)), vvec0(base_ind(swp, end), swp), '<');
                                plot(tvec0(last_ind(swp, 1)), vvec0(last_ind(swp, 1), swp), '>');
                                plot(tvec0(last_ind(swp, end)), vvec0(last_ind(swp, end), swp), '<');
                            end
                        end
                        xlabel('Time (ms)')
                        ylabel('Voltage (mV)')
                        xlim(cpwin);
                        subplot(3,1,3);
                        if ~isempty(xvec_all)
                            plot(cfit, xvec_all, yvec_all); hold on;
                            xlim(cpwin - cpwin(1));
                            ylim([(dv1(sn) - 0.5) 0.5]);
                            ax = gca;
                            xlimits = get(ax, 'Xlim');
                            line(xlimits, [dv1(sn) dv1(sn)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
                            ylimits = get(ax, 'Ylim');
                            xpos = xlimits(1) + (3/5) * (xlimits(2) - xlimits(1));
                            ypos = dv1(sn) + 1;
                            text(xpos, ypos, ['Rinput = ', num2str(Rinput(sn)), ' MOhm']);
                            xlabel('Time (ms)')
                            ylabel('Voltage (mV)')
                        end
                        figname = fullfile(outfolder, '/Rinput/', [datfn, '_Rinput', '.png']);
                        saveas(h, figname);
                        close(h);
                    end


% missing_files = {};
if ismember(datfn, broken_files) || ismember(datfn, missing_files)

%    if exist(logpath, 'file') ~= 2
%    end


                        mnoise(swp) = range(vvec3_now(hr_ind));        % range of values of the median-filtered then moving-average-filtered trace is the maximum noise

                        % Find all peaks with prominence greater than threshold (1 mV)
                        ptemp1 = find(vpeak_p > minprom);        % peak #s with prom > minprom

                        if isempty(i_first_dip_pt)
                            ioff = 0;
                        else
                            ioff = IPSC_offset(sn);
                        end
                        if ismember([datfn, '_', num2str(swp)], LTSb4spike)
                            ptemp3 = ptemp1;
                        else
                            ptemp3 = setdiff(ptemp1, ptemp2);        % peak #s that could possibly be LTS
                        end
                        % Find all peaks with the first spike occuring after the peak
                        sp2pk_i = round(sp2pk_t/sims);
                        ptemp2 = find(vpeak_i < sp1sti + sp2pk_i);    % peak #s with peak index < spike index - 10 ms
                                                % accounts for peaks with no spikes (where sp1sti == 0)

                        % Find the narrowest voltage peak that could possibly be LTS
                        [np2der(swp), pp3] = min(vpeak_2der(ptemp3));
                        psel = ptemp3(pp3);
                        % It's an LTS only if the second derivative reaches (<=) threshold
                        if np2der(swp) > lts_thr
                            ltsv(swp) = NaN;
                            ltst(swp) = NaN;
                        end

                        if ismember([datfn, '_', num2str(swp)], LTSb4spike)
                            ptemp3 = ptemp1;
                        else
                            ptemp3 = setdiff(ptemp1, ptemp2);        % peak #s that could possibly be LTS
                        end
                        % Find all peaks with the first spike occuring after the peak
                        ptemp2 = find(vpeak_i < sp1sti + sp2pk_i);    
                        end

                    fprintf('Analyzing input resistance ...\n');

%                    ft2 = fittype('-a*exp(-(x-t0)/b)-c*exp(-(x-t0)/d)+e', 'problem', 't0', 'e');    % double exponential
                            dv0(swp) = basev - lastv; % change in membrane potential on voltage trace (mV)
                            xvecs{swp} = tvec0(cpstart:cpend);
                            yvecs{swp} = vvec0(cpstart:cpend, swp) - ;
                            if dv0(swp) <= 0
                                Rinput(swp) = NaN;    % To be detected in .csv file
                            else
                                % Fit voltage trace with two exponentials
                                cfit{swp} = fit(xvecs{swp}, yvecs{swp}, ft, 'problem', tvec0(cpstart), ...
                                    'StartPoint', [dv0(swp), 2, dv0(swp), 2, basev - 2 * dv0(swp)], ...
                                    'Lower', [0, 0, 0, 0, basev - 20 * dv0(swp)], ...
                                    'Upper', [10 * dv0(swp), 100, 10 * dv0(swp), 100, basev]); 
                                    % typical tau ~ 2 ms
                                dv1(swp) = -(cfit{swp}.a + cfit{swp}.c);
                                Rinput(swp) = (dv1(swp)*10^-3)/(cpa(swp)*10^-12)/10^6;    
                                    % input resistance in MOhm
                            end
                    subplot(2,2,4);
                    hist(Rinput);
                    xlabel('Input resistance (MOhm)')
                    ylabel('# of sweeps')

                            % Find the voltage trace corresponding to the return from the current pulse
                            after_rise_pt = find(ivec1_part2(:, swp) > cpa(swp)/4, 1);
                            rtstart = (ivec1_part2_begin - 1) + after_rise_pt;    % index of first return
                            rtend = round(IPSC_start_time/sims);            % index right before IPSC
                            xvecs2{swp} = tvec0(rtstart:rtend);
                            yvecs2{swp} = vvec0(rtstart:rtend, swp);
                            cfit2{swp} = fit(xvecs2{swp}, yvecs2{swp}, ft2, 'problem', tvec0(rtstart), basev ...
                                'StartPoint', [dv0(swp), 2, dv0(swp), 2]); 
                                % typical tau ~ 2 ms

%                            Rinput(swp) = ??;    % input resistance in MOhm

                            ddv3_now_rest = ddv3_now(ind4:end);    % the rest of ddv3_now
                            ind5 = find(ddv3_now_rest > lts_thr, 1);    % first point that returns to above threshold
                            [np2der(swp) ind6] = min(ddv3_now(ind4:((ind4 - 1) + ind5)));
                            itemp2 = (ind4 - 1) + ind6;

                            ind7 = find(vvec0_now(1:npi(swp)) < bl_mean(swp), 1, 'last');
                                                % last point below baseline
                            if ind7 == npi(swp)            % entire LTS below baseline, no spikes possible
                                btime(swp) = NaN;
                                spb(swp) = NaN;
                            else
                                ind8 = find(vvec0_now(npi(swp):end) < bl_mean(swp), 1);
                                                    % first point below baseline
                                bon_lbi(swp) = ind7;
                                bon_ubi(swp) = (npi(swp) - 1) + ind8;
                                    ind9 = find(vvec0_now(bon_lbi(swp):allspi{swp}(1)) > sp_thr, 1);
                                                    % first point greater than spike threshold
                            end

                            ind4 = find(v_now(1:vpeak_i(p)) < vpeak_a(p) - vpeak_p(p), 1, 'last');
                            if isempty(ind4)
                                vpeak_lb(p) = 1;
                            else
                                vpeak_lb(p) = ind4;
                            end
                            ind5 = find(v_now(vpeak_i(p):end) < vpeak_a(p) - vpeak_p(p), 1);
                            if isempty(ind5)
                                vpeak_ub(p) = length(v_now);
                            else
                                vpeak_ub(p) = ind5;
                            end

                                ind6 = find(vvec0_now(bon_lbi(swp):allspi{swp}(1)) ...
                                         < bpeaks_a(1) - bpeaks_p(1), 1, 'last');


fprintf('Finding IPSC peaks ...\n');
IPSC_ind = round(ipscpwin(1)/sims);        % Assume no IPSC offset
ind = IPSC_ind:round(ipscpwin(2)/sims);        % indices of interest
ivec0_part1 = ivec0(ind, :);            % Use original current trace
ivec0_part1_begin = ind(1);
ipeak = zeros(1, nswps);
ipeak_ind = zeros(1, nswps);
ipeaktime = zeros(1, nswps);
parfor swp = 1:nswps                % FOR each sweep
    [ipeak(swp), itemp1] = min(ivec0_part1(:, swp));
    ipeak_ind(swp) = (ivec0_part1_begin - 1) + itemp1;
    ipeaktime(swp) = (ipeak_ind(swp) - IPSC_ind) * sims;
end
if plotIPSCpeakflag == 1
    h = figure('Visible', 'off');
    set(h, 'Name', 'IPSC peak amplitude analysis');
    clf(h);
    for swp = 1:nswps            % Plot each current trace and mark peak amplitude
        plot(tvec0(ind), ivec0(ind, swp)); hold on; 
        plot(tvec0(ipeak_ind(swp)), ipeak(swp), 'Marker', 'x', 'MarkerSize', 12);
    end
    xlabel('Time (ms)')
    ylabel('Current (pA)')
    title(['IPSC peak amplitude analysis for ', datfn]);
    figname = fullfile(outfolder, '/IPSCpeak/', [datfn, '_IPSCpeak', '.png']);
    saveas(h, figname);
    close(h);
end


% obsolete
minprom = 1;    % minimum LTS prominence in mV
maxnoiseprom1 = 0;    % maximum noise prominence for detecting narrowest peak in mV
maxnoiseprom2 = 2;    % maximum noise prominence for updating peak boundaries in mV


blw = 20;    % width in ms for calculating baseline voltage (holding potential)
lts_thr = -0.0023;        % 2nd derivative in V^2/s^2 below which defines an LTS peak
lts_thr_alt = -0.0086823;    % 2nd derivative in V^2/s^2 above which is the "gray area"
sp_thr = -30;    % amplitude threshold in mV for detecting a spike (the highest LTS peak is -34.01 mV)
sp2pk_t = 0;    % minimum time from the first spike to the peak of the LTS (ms)


fprintf('Finding LTSs ...\n');
hr_ind = round(hrange(1)/sims):1:round(hrange(2)/sims);
                    % indices for calculating maxnoise
ndp_mafw2 = round(mafw2/sims);
bl_ind = IPSC_ind - fliplr(1:round(blw/sims));    
                    % indices for calculating baseline voltage
dvvec3 = zeros(ndps - 1, nswps);    % derivative of vvec3
dvvec3_sm = zeros(ndps - 1, nswps);    % moving-average-filtered dvvec3
ddvvec3 = zeros(ndps - 2, nswps);    % derivative of dvvec3_sm
ind3 = cell(1, nswps);            % indices of interest
npi = zeros(1, nswps);            % narrowest peak index
np_lbi_old = zeros(1, nswps);        % narrowest peak lower bound index
np_ubi_old = zeros(1, nswps);        % narrowest peak upper bound index
np_lbi = zeros(1, nswps);        % narrowest peak lower bound index, updated
np_ubi = zeros(1, nswps);        % narrowest peak upper bound index, updated
isspontaneous = zeros(1, nswps);    % 1 if it's a spontaneous spike
mxsli = zeros(1, nswps);        % LTS maximum slope index
bon_i = zeros(1, nswps);        % Burst onset index
allspi = cell(1, nswps);        % All spike indices (could be burst or spontaneous spike)
pk_class_label = cell(1, nswps);    % Classification notes for selected peak
ind3_end = round(ltswin(2)/sims);

parfor swp = 1:nswps
    % Calculate maxnoise
    mnoise(swp) = range(vvec3_now(hr_ind));        % range of values of the median-filtered then moving-average-filtered trace is the maximum noise

    % Calculate baseline voltage (holding potential)
    bl_mean(swp) = mean(vvec1_now(bl_ind));        % baseline voltage (actual holding potential)
                            % Previously uses vvec0_now

    % Set up 2nd derivative of median-filtered voltage trace
    dvvec3_now = diff(vvec3_now)./diff(tvec0);    % differentiate vvec3, the median-filtered then smoothed voltage trace
    dvvec3(:, swp) = dvvec3_now;            % store dvvec3 for plotting
    dvvec3_sm_now = smooth(dvvec3_now, ndp_mafw2);    % smooth out dvvec3
    dvvec3_sm(:, swp) = dvvec3_sm_now;        % store dvvec3_sm for plotting
    ddvvec3_now = diff(dvvec3_sm_now)./diff(tvec0(2:end));    % differentiate dvvec3_sm
    ddvvec3(:, swp) = ddvvec3_now;            % store ddvvec3 for plotting
    ind3{swp} = ipeak_ind(swp):ind3_end;
    v3_now = vvec3_now(ind3{swp});            % voltage vector of interest for detecting LTS
    v0_now = vvec0_now(ind3{swp});            % voltage vector of interest for detecting spikes (ap)
    v_now_begin = ipeak_ind(swp);
    ddv3_now = ddvvec3_now(ind3{swp} - 1);        % 2nd derivative vector of interest
                            % two differentiations: 
                            % shift right (don't -1) than shift left (-1)

    % Find all voltage peaks, locate peak bounds, compute most negative 2nd derivatives, detect spikes
    [vpeak_a, vpeak_i, vpeak_w, vpeak_p] = findpeaks(v3_now);    % find all voltage peaks
    npks = length(vpeak_a);
    vpeak_lb = zeros(npks, 1);
    vpeak_ub = zeros(npks, 1);
    vpeak_2der = zeros(npks, 1);
    v0_pk = cell(npks, 1);
    v0_pk_begin = zeros(npks, 1);
    spi = cell(npks, 1);
    sp1sti = zeros(npks, 1);
    for p = 1:npks
        % Find peak lower bounds
        if vpeak_i(p) < 3            % findpeaks will not work for < 3 points
            vpeak_lb(p) = 1;
        else
            [amp4, ind4] = findpeaks(-flipud(v3_now(1:vpeak_i(p))), ...
                    'MinPeakProminence', 0, 'NPeaks', 1);
                % first minimum to the left: flip and invert, 
                % then find first voltage peak
            if isempty(ind4)
                vpeak_lb(p) = 1;
            else
                vpeak_lb(p) = (vpeak_i(p) + 1) - ind4;
            end
        end
        % Find peak upper bounds
        if vpeak_i(p) > length(v3_now) - 2    % findpeaks will not work for < 3 points
            vpeak_ub(p) = length(v3_now);
        else
            [amp5, ind5] = findpeaks(-v3_now(vpeak_i(p):end), ...
                    'MinPeakProminence', 0, 'NPeaks', 1);
                % first minimum to the right: invert, 
                % then find first voltage peak
            if isempty(ind5)
                vpeak_ub(p) = length(v3_now);
            else
                vpeak_ub(p) = (vpeak_i(p) - 1) + ind5;
            end
        end
        % Find most negative 2nd derivative over the entire peak
        if vpeak_lb(p) == vpeak_ub(p)
            vpeak_2der(p) = ddv3_now(vpeak_lb(p));
        else
            vpeak_2der(p) = min(ddv3_now(vpeak_lb(p):vpeak_ub(p)));
        end
        % Detect spikes in original trace within the peak
        v0_pk{p} = v0_now(vpeak_lb(p):vpeak_ub(p));        % take only the peak part
        v0_pk_begin(p) = vpeak_lb(p);
        [pspikes_a, pspikes_i] = findpeaks(v0_pk{p});    % find all "spikes" within the peak
        stemp1 = find(pspikes_a > sp_thr);    % action potentials must be greater than threshold
        if ~isempty(stemp1)
            % Record spike indices relative to v0_now or v3_now
            spi{p} = (v0_pk_begin(p) - 1) + pspikes_i(stemp1);
            % Record index of first spike relative to v0_now or v3_now
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
    ptemp1 = find(vpeak_p > mnoise(swp));    % peak #s with prom > maxnoise         
                % (previously 1 mV, now different for each trace)

    % Find all peaks with the second derivative reaching (<=) threshold
    ptemp2 = find(vpeak_2der <= lts_thr);

    % Find all peaks that are LTS candidates by threshold
    ptemp3 = intersect(ptemp1, ptemp2);    % peak #s that are LTSs by threshold

    if isempty(ptemp1)        % Condition (3)
        % Not LTS by prominence
        ltsv(swp) = NaN;
        ltst(swp) = NaN;
        mxslt(swp) = NaN;
        mxslv(swp) = NaN;

        % Find the narrowest voltage peak
        [np2der(swp), psel] = min(vpeak_2der);    % find the narrowest voltage peak
        pk_class(swp) = 1;
        pk_class_label{swp} = ['not LTS: peak prominence ', ...
                    num2str(vpeak_p(psel)), ...
                    ' mV <= maximum noise ', ...
                    num2str(mnoise(swp)), ' mV'];
    elseif isempty(ptemp3)        % Condition (2)
        % Not LTS by narrowness
        ltsv(swp) = NaN;
        ltst(swp) = NaN;
        mxslt(swp) = NaN;
        mxslv(swp) = NaN;

        % Find the narrowest voltage peak with prominence greater than maximum noise
        [np2der(swp), pp1] = min(vpeak_2der(ptemp1));
        psel = ptemp1(pp1);
        pk_class(swp) = 2;
        pk_class_label{swp} = ['not LTS: 2nd derivative ', ...
                    num2str(np2der(swp)), ...
                    ' V^2/s^2 > LTS threshold ', ...
                    num2str(lts_thr), ' V^2/s^2'];
    else                % Condition (1)
        % Select the first peak that is an LTS candidate by prominence & 2nd der
        psel = ptemp3(1);
        np2der(swp) = vpeak_2der(psel);

        % Check whether it's a spontaneous spike
        %     (Based on following observation of shape:
        %         LTS:         first spike occurs before "LTS" peak on mfmaf trace, 
        %            except in two cases (A092910_0001_23 & F101210_0000_1),
        %             where the first spike occurred after the peak
        %     spontaneous spike:     first spike occurs after "LTS" peak on mfmaf trace)
        sp2pk_i = round(sp2pk_t/sims);    % obsolete: currently set to 0
        if ~ismember(filebase, LTSb4spike) ...
            && sp1sti(psel) ~= 0 ...
            && vpeak_i(psel) < sp1sti(psel) + sp2pk_i
            % Not LTS by shape
            ltsv(swp) = NaN;
            ltst(swp) = NaN;
            mxslt(swp) = NaN;
            mxslv(swp) = NaN;

            isspontaneous(swp) = 1;
            pk_class(swp) = 3;
            pk_class_label{swp} = ['not LTS: peak index ', ...
                        num2str(vpeak_i(psel)), ...
                        ' < index of first spike ', ...
                        num2str(sp1sti(psel))];
        end
    end

    % Record prominence of selected peak
    pk_prom(swp) = vpeak_p(psel);

    % Record indices relative to tvec0 or vvec0_now
    npi(swp) = (v_now_begin - 1) + vpeak_i(psel);        % narrowest peak index
    npt(swp) = (npi(swp) - IPSC_ind) * sims;        % narrowest peak time (delay)
    % The following may be changed later for bursts
    np_lbi(swp) = (v_now_begin - 1) + vpeak_lb(psel);    % narrowest peak lower bound index
    np_ubi(swp) = (v_now_begin - 1) + vpeak_ub(psel);    % narrowest peak upper bound index

    % Save old peak boundary indices for plotting
    np_lbi_old(swp) = np_lbi(swp);
    np_ubi_old(swp) = np_ubi(swp);

    % Record spike indices relative to vvec0_now & count spikes per peak
    if ~isempty(spi{psel})
        allspi{swp} = (v_now_begin - 1) + spi{psel};    % spike indices
        spp(swp) = length(spi{psel});            % count spikes per peak
    end

    % Find LTS amplitude & delay
    %                        if ~isempty(ptemp3) && np2der(swp) <= lts_thr
    if ~isempty(ptemp3) && isspontaneous(swp) ~= 1
        ltsv(swp) = vvec1_now(npi(swp));    % LTS amplitude, use median-filtered voltage trace
        ltst(swp) = npt(swp);            % LTS delay

        [mxslv(swp), ind1] = max(dvvec3_sm_now(np_lbi(swp):np_ubi(swp)));    
                            % maximum slope in V/s
        mxsli(swp) = np_lbi(swp) + ind1 - 1;
                            % index of maximum slope
        mxslt(swp) = (mxsli(swp) - IPSC_ind) * sims;
                            % time in ms of maximum slope after IPSC starts

    end

    % Determine whether it's a "burst" and count action potentials (spikes)
    if isnan(ltst(swp))                % must be an "LTS" to begin with
        btime(swp) = NaN;
        spb(swp) = NaN;
    else
        if isempty(allspi{swp})            % no spikes detected, not a "burst"
            btime(swp) = NaN;
            spb(swp) = NaN;
            if np2der(swp) <= lts_thr_alt
                pk_class(swp) = 6;
                pk_class_label{swp} = ['LTS with no burst; ', ...
                            'definite: 2nd der ', ...
                            num2str(np2der(swp)), ...
                            ' V^2/s^2'];
            else
                pk_class(swp) = 4;
                pk_class_label{swp} = ['LTS with no burst; ', ...
                            'contentious: 2nd der ', ...
                            num2str(np2der(swp)), ...
                            ' V^2/s^2 > ', ...
                            'alt LTS thr ', ...
                            num2str(lts_thr_alt), ...
                            ' V^2/s^2'];
            end
        else
            if np2der(swp) <= lts_thr_alt
                pk_class(swp) = 7;
                pk_class_label{swp} = ['LTS with burst; definite: ', ...
                            '2nd der ', ...
                            num2str(np2der(swp)), ...
                            ' V^2/s^2'];
            else
                pk_class(swp) = 5;
                pk_class_label{swp} = ['LTS with burst; contentious: ', ...
                            '2nd der ', ...
                            num2str(np2der(swp)), ...
                            ' V^2/s^2 > ', ...
                            'alt LTS thr ', ...
                            num2str(lts_thr_alt), ...
                            ' V^2/s^2'];
            end

            % Re-detect spikes by re-finding peak bounds using MinPeakProminence
            % Update peak lower bounds
            if vpeak_i(psel) < 3            % findpeaks will not work for < 3 points
                vpeak_lb(psel) = 1;
            else
                [amp4, ind4] = findpeaks(-flipud(v3_now(1:vpeak_i(psel))), ...
                        'MinPeakProminence', mnoise(swp), 'NPeaks', 1);
                    % first minimum to the left: flip and invert, 
                    % then find first voltage peak with prominence >= maxnoise (previously 2 mV, now different for each trace)
                if isempty(ind4)
                    vpeak_lb(psel) = 1;
                else
                    vpeak_lb(psel) = (vpeak_i(psel) + 1) - ind4;
                end
            end
            % Update peak upper bounds
            if vpeak_i(psel) > length(v3_now) - 2    % findpeaks will not work for < 3 points
                vpeak_ub(psel) = length(v3_now);
            else
                [amp5, ind5] = findpeaks(-v3_now(vpeak_i(psel):end), ...
                        'MinPeakProminence', mnoise(swp), 'NPeaks', 1);
                    % first minimum to the right: invert, 
                    % then find first voltage peak with prominence >= maxnoise (previously 2 mV, now different for each trace)
                if isempty(ind5)
                    vpeak_ub(psel) = length(v3_now);
                else
                    vpeak_ub(psel) = (vpeak_i(psel) - 1) + ind5;
                end
            end
            % Detect spikes in original trace within the peak
            v0_pk{psel} = v0_now(vpeak_lb(psel):vpeak_ub(psel));        % take only the peak part
            v0_pk_begin(psel) = vpeak_lb(psel);
            [pspikes_a, pspikes_i] = findpeaks(v0_pk{psel});    % find all "spikes" within the peak
            stemp1 = find(pspikes_a > sp_thr);    % action potentials must be greater than threshold
            if ~isempty(stemp1)
                % Record spike indices relative to v0_now or v3_now
                spi{psel} = (v0_pk_begin(psel) - 1) + pspikes_i(stemp1);
            end

            % Update indices relative to tvec0 or vvec0_now
            npi(swp) = (v_now_begin - 1) + vpeak_i(psel);        % narrowest peak index
            npt(swp) = (npi(swp) - IPSC_ind) * sims;        % narrowest peak time (delay)
            np_lbi(swp) = (v_now_begin - 1) + vpeak_lb(psel);    % narrowest peak lower bound index
            np_ubi(swp) = (v_now_begin - 1) + vpeak_ub(psel);    % narrowest peak upper bound index

            % Update spike indices and spikes per peak
            if ~isempty(spi{psel})
                allspi{swp} = (v_now_begin - 1) + spi{psel};    % spike indices
                spp(swp) = length(spi{psel});            % count spikes per peak
            end

            % Spikes per burst is the the same spikes per peak but not zero
            spb(swp) = spp(swp);        % spikes per burst

            % Find burst onset time (delay)
            vvec0_b = v0_pk{psel};        % voltage trace of burst
            bspk1i = sp1sti(psel) - v0_pk_begin(psel) + 1;    % index of first spike in vvec0_b

            % Find first minimum to the left of first spike
            if bspk1i < 3    % findpeaks will not work for < 3 points
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
            if bspk1i > length(vvec0_b) - 2    % findpeaks will not work for < 3 points
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
            amp8 = max([-amp6 -amp7]);    % take the higher of the two minimums

            % Burst onset index is the last point
            %     lower than the base of the first spike
            %     Find index in terms of v0_now
            ind8 = find(v0_now(1:sp1sti(psel)) < amp8, 1, 'last');
            % Convert index in terms of vvec0_now
            if isempty(ind8)
                bon_i(swp) = (v_now_begin - 1) + v0_pk_begin(psel);
            else
                bon_i(swp) = (v_now_begin - 1) + ind8;
            end
            btime(swp) = (bon_i(swp) - IPSC_ind)*sims;
                            % burst onset time (delay)
        end
    end
end

if plotLTSflag == 1
    fprintf('Plotting voltage traces ...\n');
    bl_start = bl_ind(1);
    bl_end = bl_ind(end);
    thiscellname = cellnames{celln};
    parfor swp = 1:nswps                % FOR each sweep
        vvec0_now = vvec0(:, swp);
        vvec1_now = vvec1(:, swp);
        vvec2_now = vvec2(:, swp);
        vvec3_now = vvec3(:, swp);
        dvvec3_now = dvvec3(:, swp);
        dvvec3_sm_now = dvvec3_sm(:, swp);
        ddvvec3_now = ddvvec3(:, swp);
        filebase = [datfn, '_', num2str(swp)];

        % Plot voltage traces
        h = figure(swp*10 + 1);
%                        h = figure('Visible', 'off');
        set(h, 'Name', 'Voltage traces');
        clf(h);
        plot(tvec0, vvec0_now, 'b-', 'LineWidth', 0.5); hold on
        plot(tvec0, vvec1_now, 'g-', 'LineWidth', 0.5);
        plot(tvec0, vvec3_now, 'r-', 'LineWidth', 0.5);
        if ~isnan(ltst(swp)) && ~isnan(btime(swp))    % LTS with bursts
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'go', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'gd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
            plot(tvec0(bon_i(swp)), vvec0_now(bon_i(swp)), 'g>', 'MarkerSize', 10);
            plot(tvec0(allspi{swp}), vvec0_now(allspi{swp}), 'gx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                pk_class_label{swp}, 'burst onset', 'spikes', ...
                'Location', 'SouthOutside')
        elseif ~isnan(ltst(swp))            % LTS without bursts
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'bo', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'bd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
            legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                pk_class_label{swp}, 'Location', 'SouthOutside')
        else                        % not LTS
            if np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'rx', 'MarkerSize', 10);
            else
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'kx', 'MarkerSize', 10);
            end
            if isempty(allspi{swp})            % noise
                legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                    pk_class_label{swp}, 'Location', 'SouthOutside')
            else                    % spontaneous spikes
                plot(tvec0(allspi{swp}(1)), vvec0_now(allspi{swp}(1)), 'rx', 'MarkerSize', 10);
                legend('raw trace', 'median-filtered', 'median-filtered then moving-average-filtered', ...
                    pk_class_label{swp}, 'spontaneous spike', ...
                    'Location', 'SouthOutside')
            end
        end
        plot(tvec0(bl_start), vvec1_now(bl_start), 'g>');
        plot(tvec0(bl_end), vvec1_now(bl_end), 'y<');
        plot(tvec0(np_lbi_old(swp)), vvec3_now(np_lbi_old(swp)), 'k*');
        plot(tvec0(np_ubi_old(swp)), vvec3_now(np_ubi_old(swp)), 'k*');
        plot(tvec0(np_lbi(swp)), vvec3_now(np_lbi(swp)), 'r*');
        plot(tvec0(np_ubi(swp)), vvec3_now(np_ubi(swp)), 'r*');
        xlim([tvec0(1) tvec0(end)]);
        xlabel('Time (ms)');
        ylabel('Voltage (mV)');
        title(strrep(filebase, '_', '\_'));
        figname = fullfile(outfolder, '/vtraces/', [filebase, '.png']);
        saveas(h, figname);
        figure(h);
        ylim([-120 20]);
        figname = fullfile(outfolder, '/vtraces_scaled/', [filebase, '_scaled.png']);
        saveas(h, figname);
        if np2der(swp) > lts_thr_alt ...
            && np2der(swp) <= lts_thr    % in "gray area"
            figname = fullfile(outfolder, '/gray_area_traces/', [filebase, '_scaled.png']);
            saveas(h, figname);
        end
%                        close(h);

        % Plot LTS analysis
        xlimits = ltswin;
        h = figure(swp*10 + 2);
%                         h = figure('Visible', 'off');
        set(h, 'Name', 'LTS analysis, moving-average-filtered trace');
        clf(h);
        subplot(3,1,1) % voltage trace
        plot(tvec0, vvec3_now, 'r-', 'LineWidth', 0.5); hold on
        if ~isnan(ltst(swp)) && ~isnan(btime(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'go', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'gd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
        elseif ~isnan(ltst(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'bo', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'bd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
        else
            if np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'rx', 'MarkerSize', 10);
            else
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'kx', 'MarkerSize', 10);
            end
        end
        plot(tvec0(np_lbi(swp)), vvec3_now(np_lbi(swp)), 'r*');
        plot(tvec0(np_ubi(swp)), vvec3_now(np_ubi(swp)), 'r*');
        plot(tvec0(np_lbi_old(swp)), vvec3_now(np_lbi_old(swp)), 'k*');
        plot(tvec0(np_ubi_old(swp)), vvec3_now(np_ubi_old(swp)), 'k*');
        xlim(xlimits);
        xlabel('Time (ms)');
        ylabel('Voltage (mV)');
        title(['LTS analysis for ', strrep(filebase, '_', '\_'), ', moving-average-filtered trace']);
        subplot(3,1,2) % 1st derivative of voltage trace
        plot(tvec0(2:end), dvvec3_now, 'k-', 'LineWidth', 0.5); hold on
        plot(tvec0(2:end), dvvec3_sm_now, 'r-', 'LineWidth', 0.5);
        if ~isnan(ltst(swp)) && ~isnan(btime(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), dvvec3_now(npi(swp)), 'go', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), dvvec3_now(mxsli(swp)), 'gd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), dvvec3_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), dvvec3_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
        elseif ~isnan(ltst(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), dvvec3_now(npi(swp)), 'bo', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), dvvec3_now(mxsli(swp)), 'bd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), dvvec3_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), dvvec3_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
        else
            if np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), dvvec3_now(npi(swp)), 'rx', 'MarkerSize', 10);
            else
                plot(tvec0(npi(swp)), dvvec3_now(npi(swp)), 'kx', 'MarkerSize', 10);
            end
        end
        plot(tvec0(np_lbi(swp)), dvvec3_now(np_lbi(swp)), 'r*');
        plot(tvec0(np_ubi(swp)), dvvec3_now(np_ubi(swp)), 'r*');
        plot(tvec0(np_lbi_old(swp)), dvvec3_now(np_lbi_old(swp)), 'k*');
        plot(tvec0(np_ubi_old(swp)), dvvec3_now(np_ubi_old(swp)), 'k*');
        xlim(xlimits);
        xlabel('Time (ms)');
        ylabel('dV/dT');
        legend('unsmoothed', 'smoothed');
        subplot(3,1,3) % 2nd derivative of voltage trace
        plot(tvec0(3:end), ddvvec3_now, 'k-', 'LineWidth', 0.5); hold on
        line(xlimits, [lts_thr lts_thr], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
        line(xlimits, [lts_thr_alt lts_thr_alt], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 0.5);
        xlim(xlimits);
        ax = gca;
        ylimits = get(ax, 'YLim');
        line([tvec0(ipeak_ind(swp)) tvec0(ipeak_ind(swp))], ...
            ylimits, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
        if ~isnan(ltst(swp)) && ~isnan(btime(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), ddvvec3_now(npi(swp)), 'go', 'MarkerSize', 10);
            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), ddvvec3_now(npi(swp)), 'ro', 'MarkerSize', 10);
            end
        elseif ~isnan(ltst(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), ddvvec3_now(npi(swp)), 'bo', 'MarkerSize', 10);
            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), ddvvec3_now(npi(swp)), 'ro', 'MarkerSize', 10);
            end
        else
            if np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), ddvvec3_now(npi(swp)), 'rx', 'MarkerSize', 10);
            else
                plot(tvec0(npi(swp)), ddvvec3_now(npi(swp)), 'kx', 'MarkerSize', 10);
            end
        end
        plot(tvec0(np_lbi(swp)), ddvvec3_now(np_lbi(swp)), 'r*');
        plot(tvec0(np_ubi(swp)), ddvvec3_now(np_ubi(swp)), 'r*');
        plot(tvec0(np_lbi_old(swp)), ddvvec3_now(np_lbi_old(swp)), 'k*');
        plot(tvec0(np_ubi_old(swp)), ddvvec3_now(np_ubi_old(swp)), 'k*');
        xlabel('Time (ms)');
        ylabel('d2V/dT2');
        figname = fullfile(outfolder, '/LTSanalysis/', [filebase, '_LTSanalysis', '.png']);
        saveas(h, figname);
%                        close(h);

        % Plot burst analysis
        h = figure(swp*10 + 3);
%                        h = figure('Visible', 'off');
        set(h, 'Name', 'Burst analysis, original trace');
        clf(h);
        plot(tvec0(np_lbi(swp):np_ubi(swp)), vvec0_now(np_lbi(swp):np_ubi(swp)), 'b-', 'LineWidth', 0.5); hold on
        plot(tvec2, vvec2_now, 'g.', 'LineWidth', 0.5);
        plot(tvec0(np_lbi(swp):np_ubi(swp)), vvec3_now(np_lbi(swp):np_ubi(swp)), 'r-', 'LineWidth', 0.5);
        if ~isnan(ltst(swp)) && ~isnan(btime(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'go', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'gd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
            plot(tvec0(bon_i(swp)), vvec0_now(bon_i(swp)), 'g>', 'MarkerSize', 10);
            plot(tvec0(allspi{swp}), vvec0_now(allspi{swp}), 'gx', 'MarkerSize', 10);
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                pk_class_label{swp}, 'burst onset', 'spikes', ...
                'Location', 'SouthOutside')
        elseif ~isnan(ltst(swp))
            if np2der(swp) <= lts_thr_alt
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'bo', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'bd', 'MarkerSize', 8);

            elseif np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'ro', 'MarkerSize', 10);
                plot(tvec0(mxsli(swp)), vvec1_now(mxsli(swp)), 'rd', 'MarkerSize', 8);

            end
            legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                pk_class_label{swp}, 'Location', 'SouthOutside')
        else
            if np2der(swp) > lts_thr_alt ...
                && np2der(swp) <= lts_thr    % in "gray area"
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'rx', 'MarkerSize', 10);
            else
                plot(tvec0(npi(swp)), vvec1_now(npi(swp)), 'kx', 'MarkerSize', 10);
            end
            if isempty(allspi{swp})
                legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                    pk_class_label{swp}, 'Location', 'SouthOutside')
            else
                plot(tvec0(allspi{swp}(1)), vvec0_now(allspi{swp}(1)), 'rx', 'MarkerSize', 10);
                legend('raw trace', 'median-filtered then resampled', 'median-filtered then moving-average-filtered', ...
                    pk_class_label{swp}, 'spontaneous spike', ...
                    'Location', 'SouthOutside')
            end
        end
        plot(tvec0(np_lbi_old(swp)), vvec3_now(np_lbi_old(swp)), 'k*');
        plot(tvec0(np_ubi_old(swp)), vvec3_now(np_ubi_old(swp)), 'k*');
        plot(tvec0(np_lbi(swp)), vvec3_now(np_lbi(swp)), 'r*');
        plot(tvec0(np_ubi(swp)), vvec3_now(np_ubi(swp)), 'r*');
        ax = gca;
        if isspontaneous(swp) == 1
            text((19/20)*ax.XLim(1) + (1/20)*ax.XLim(2), ...
                (1/20)*ax.YLim(1) + (19/20)*ax.YLim(2), ...
                'Spontaneous spike!');
        end

        xlim([np_lbi(swp) np_ubi(swp)] * sims);
        xlabel('Time (ms)');
        ylabel('Voltage (mV)');
        title(['Burst analysis for ', strrep(filebase, '_', '\_'), ', original trace'])
        figname = fullfile(outfolder, '/burstanalysis/', [filebase, '_burstanalysis', '.png']);
        saveas(h, figname);
%                        close(h);

    end
end

 || ismember(datfn, missing_files)

                            % Record current and voltage traces; 
                            % NOTE: some data points seem corrupted
                            if ismember(datfn, flipped_files)
                                gvec0(:, swp) = d(:, 2, swp);    % conductance trace in pS OR nS
                                ivec0(:, swp) = d(:, 3, swp);    % current trace in pA
                            else
                                gvec0(:, swp) = d(:, 3, swp);    % conductance trace in pS OR nS
                                ivec0(:, swp) = d(:, 2, swp);    % current trace in pA
                            end

                            % Convert all conductance traces to nS
                            if max(gvec0(:, swp)) > 100            % was in pS
                                gvec0(:, swp) = gvec0(:, swp) / 1000;    % now in nS
                            end

                            % Fix current and conductance traces that were arbitrarily scaled
                            cpmidind = round(cpmid/sims):round((cpmid + mvw)/sims);
                                            % indices for measuring cp peak
                            ivec0_now = ivec0(:, swp);
                            cpa_ap(swp) = mean(ivec0_now(cpmidind));    % approximate cp amplitude
                            if cpa_ap(swp) < -105 && cpa_ap(swp) > -145    % about -125 pA
                                ivec0(:, swp) = ivec0(:, swp) / 2.5; % arbitrarily scaled by 2.5
                                gvec0(:, swp) = gvec0(:, swp) / 2.5; % arbitrarily scaled by 2.5
                            elseif cpa_ap(swp) < -230 && cpa_ap(swp) > -270  % about -250 pA
                                ivec0(:, swp) = ivec0(:, swp) / 5; % arbitrarily scaled by 5
                                gvec0(:, swp) = gvec0(:, swp) / 5; % arbitrarily scaled by 5
                            end


                            % Convert all conductance traces to nS
                            %%% 2016-11-07 THIS DOESN'T CORRECT ALL FILES!!
                            if max(gvec0_temp) > 100            % was in pS
                                gvec0_temp = gvec0_temp / 1000;        % now in nS
                            end


row = 0;        % Current row # minus one in the .csv file (each sweep gets a row); first row is for the log header

        % Used only if appending new data
        rowthisp = row;

        % Print properties of each sweep into the datalog
        % Must be consistent with logheader & logvariables
        if saveswpinfoflag == 1
            fid = fopen(logpath, 'a');
            for k = (rowthisp+1):row
                fprintf(fid, ['%s, %d, ' ...
                '%d, %d, %d, %d, ' ...
                '%g, %g, %g, ', ...
                '%g, %g, ' ...
                '%10.5f, %10.5f, %10.5f, %10.5f, %10.5f, ' ...
                '%10.5f, %10.5f, %10.5f, %10.5f, ' ...
                '%10.5f, %10.5f, %d, %d, ' ...
                '%10.5f, %10.5f, %10.5f, %10.5f, ' ...
                '%10.5f, %d, %10.5f, %10.5f', ...
                '%10.5f, %10.5f\n'], ...
                fnrow{k}, cellidrow(k), ...
                prow(k), vrow(k), grow(k), swpnrow(k), ...
                gabab_amp(k), gabab_Trise(k), gabab_TfallFast(k), ...
                gabab_TfallSlow(k), gabab_w(k), ...
                currpulse(k), Rin(k), ioffset_old(k), imint(k), imin(k), ...
                actVhold(k), maxnoise(k), peaktime(k), peak2ndder(k), ...
                peakprom(k), peakwidth(k), peakclass(k), spikesperpeak(k), ...
                ltspeaktime(k), ltspeakval(k), maxslopetime(k), maxslopeval(k), ...
                bursttime(k), spikesperburst(k), maxspikeamp(k), minspikeamp(k), ...
                spikefrequency(k), spikeadaptation(k));
            end
            fclose(fid);
        end

        if size(d, 3) ~= length(gtag)

    save(swpdatafn, 'logheader', 'logvariables', ...
            'cellnames', 'abffullfn', 'nswps', ...
            'fnrow', 'cellidrow', ...
            'prow', 'vrow', 'grow', 'swpnrow', ...
            'gabab_amp', 'gabab_Trise', 'gabab_TfallFast', ...
            'gabab_TfallSlow', 'gabab_w', ...
            'currpulse', 'Rin', 'ioffset_old', 'imint', 'imin', ...
            'actVhold', 'maxnoise', 'peaktime', 'peak2ndder', ...
            'peakprom', 'peakwidth', 'peakclass', 'spikesperpeak', ...
            'ltspeaktime', 'ltspeakval', 'maxslopetime', 'maxslopeval', ...
            'bursttime', 'spikesperburst', 'maxspikeamp', 'minspikeamp', ...
            'spikefrequency', 'spikeadaptation', ...
            '-v7.3');


[ndps, sims, tvec0, gvec0, ivec0, vvec0, ...
    gvec1, ivec1, vvec1, ...
    ndps2, sims2, tvec2, gvec2, ivec2, vvec2, ...
    vvec3] = ...
    load_matfiles_part (newinfolder, datfn, nswps);

firstfile = fullfile(newinfolder, [datfn, '_1.mat']);
if exist(firstfile, 'file') ~= 2
    error(['This mat file: ', firstfile, ' is missing!!']);
end

parfor swp = 1:nswps        % FOR each sweep
    filebase = [datfn, '_', num2str(swp)];
    thisfile = fullfile(newinfolder, [filebase, '.mat']);
    if ~exist(thisfile, 'file')
        error(['This mat file: ', thisfile, ' is missing!!']);
    end
end

if ~isfolder(outfolder)
    mkdir(outfolder);
end
if ~isfolder(fullfile(outfolder, '/passive/'))
    mkdir(fullfile(outfolder, '/passive/'));
end
if ~isfolder(fullfile(outfolder, '/IPSCoffset_old/'))
    mkdir(fullfile(outfolder, '/IPSCoffset_old/'));
end
if ~isfolder(fullfile(outfolder, '/matfiles/'))
    mkdir(fullfile(outfolder, '/matfiles/'));
end
if ~isfolder(fullfile(outfolder, '/IPSCpeak/'))
    mkdir(fullfile(outfolder, '/IPSCpeak/'));
end
if ~isfolder(fullfile(outfolder, '/vtraces/'))
    mkdir(fullfile(outfolder, '/vtraces/'));
end
if ~isfolder(fullfile(outfolder, '/LTSanalysis/'))
    mkdir(fullfile(outfolder, '/LTSanalysis/'));
end
if ~isfolder(fullfile(outfolder, '/burstanalysis/'))
    mkdir(fullfile(outfolder, '/burstanalysis/'));
end
if ~isfolder(fullfile(outfolder, '/vtraces_scaled/'))
    mkdir(fullfile(outfolder, '/vtraces_scaled/'));
end
if ~isfolder(fullfile(outfolder, '/gray_area_traces/'))
    mkdir(fullfile(outfolder, '/gray_area_traces/'));
end
if ~isfolder(fullfile(outfolder, '/LTScouldbemissed/'))
    mkdir(fullfile(outfolder, '/LTScouldbemissed/'));
end

if isfolder('/tmp/data/m3ha/')
    homedirectory = '/tmp/data/m3ha/';
elseif isfolder('/media/adamX/m3ha/')
    homedirectory = '/media/adamX/m3ha/';
elseif isfolder('/scratch/al4ng/m3ha/')
    homedirectory = '/scratch/al4ng/m3ha/';
else
    error('Valid homedirectory does not exist!');
end

%}
