function [bl_mean, mnoise, npt, np2der, pk_prom, pk_width, pk_class, spp, ltsv, ltst, mxslv, mxslt, btime, spb, spthr, fspt, lspt, maxspi, minspi, spif, spia] = find_LTSs_many_sweeps (tvec0, vvec0, IPSC_start, ipeak_time, hrange, ltswin, plotLTSflag, outfolder, datfn, tvec2, vvec1, vvec2, vvec3)
%% Calls find_LTS.m for many voltage traces
% Usage: [bl_mean, mnoise, npt, np2der, pk_prom, pk_width, pk_class, spp, ltsv, ltst, mxslv, mxslt, btime, spb, maxspi, minspi, spif, spia] = find_LTSs_many_sweeps (tvec0, vvec0, IPSC_start, ipeak_time, hrange, ltswin, plotLTSflag, outfolder, datfn, tvec2, vvec1, vvec2, vvec3)
% Arguments:    
%
% Used by:
%        /media/adamX/m3ha/data_dclamp/dclampDataExtractor.m
% Requires:
%        /home/Matlab/Adams_Functions/find_LTS.m
%
% 2016-11-07 Moved from dclampDataExtractor.m
% 2017-12-21 Changed tabs to spaces
% 2017-12-21 Added spthr, fspt, lspt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nswps = size(vvec0, 2);

bl_mean = zeros(1, nswps);      % baseline voltage
mnoise = zeros(1, nswps);       % maximum noise in mfmaf voltage trace
npt = zeros(1, nswps);          % narrowest peak time (delay)
np2der = zeros(1, nswps);       % "narrowest peak" 2nd derivative
                                % (actually first peak that crosses threshold)
pk_prom = zeros(1, nswps);      % Prominence of selected peak
pk_width = zeros(1, nswps);     % Width of selected peak at half-prominence
pk_class = zeros(1, nswps);     % Classification number for selected peak
spp = zeros(1, nswps);          % Spikes per peak (can be 0)
ltsv = zeros(1, nswps);         % LTS amplitude
ltst = zeros(1, nswps);         % LTS delay
mxslv = zeros(1, nswps);        % LTS maximum slope amplitude
mxslt = zeros(1, nswps);        % LTS maximum slope delay
btime = zeros(1, nswps);        % Burst onset time (delay)
spb = zeros(1, nswps);          % Spikes per burst (cannot be 0)
spthr = zeros(1, nswps);        % Spike threshold (mV)
fspt = zeros(1, nswps);         % First spike time (ms)
lspt = zeros(1, nswps);         % Last spike time (ms)
maxspi = zeros(1, nswps);       % Maximum spike amplitude
minspi = zeros(1, nswps);       % Minimum spike amplitude
spif = zeros(1, nswps);         % Spike frequency (Hz)
spia = zeros(1, nswps);         % Spike adaptation (%)
parfor swp = 1:nswps            % FOR each sweep
% for swp = 1:nswps             % FOR each sweep
    vvec0_now = vvec0(:, swp);              % needed for parfor
    vvec1_now = vvec1(:, swp);              % needed for parfor
    vvec2_now = vvec2(:, swp);              % needed for parfor
    vvec3_now = vvec3(:, swp);              % needed for parfor
    filebase = [datfn, '_', num2str(swp)];  % file base for this sweep

    % Find LTS info and plot figures if plotLTSflag == 1
    [bl_mean(swp), mnoise(swp), ...
    npt(swp), np2der(swp), pk_prom(swp), pk_width(swp), ...
    pk_class(swp), spp(swp), ltst(swp), ltsv(swp), ...
    mxslt(swp), mxslv(swp), btime(swp), spb(swp), ...
    spthr(swp), fspt(swp), lspt(swp), ...
    maxspi(swp), minspi(swp), spif(swp), spia(swp)] = ...
        find_LTS (tvec0, vvec0_now, IPSC_start, ...
        ipeak_time(swp), hrange, ltswin, ...
        plotLTSflag, outfolder, filebase, ...
        tvec2, vvec1_now, vvec2_now, vvec3_now);

    close all;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%