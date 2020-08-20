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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
