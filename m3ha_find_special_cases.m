function m3ha_find_special_cases(versionNumber)
%% Looks for special cases and put traces in corresponding folder
% Usage: m3ha_find_special_cases(versionNumber)
%
% Requires:    
%       cd/create_subdir_copy_files.m

% File History:
% 2016-10-19 Created
% 2016-12-05 Moved code to create_subdir_copy_files()
% 2016-12-05 Added All_Spontaneous_Spikes, Not_LTS_by_prom && Long_Latency_LTS
% 2016-12-06 Moved code to create_subdir_copy_files.m
% 2017-03-29 - BT - Changed varargin to read version number for new vtraces folder, updated to new file hierarchy

%% Parameters used in the analyses
lts_thr = -0.0023;        % 2nd derivative in kV/s^2 below which defines an LTS peak
lts_thr_alt = -0.0081823;    % 2nd derivative in kV/s^2 above which is the "gray area"

%% Set folders for reading and saving files
inFolder = '//media/adamX/m3ha/data_dclamp/take4/backup/';
outFolder = '//media/adamX/m3ha/data_dclamp/take4/special_cases/';

%% Specify which matfile to use; assumed to be in inFolder
filetouse = ['dclampdatalog_take4_old', num2str(versionNumber), '.mat'];

%% Specify which types of figures to copy
toCopySuffix = {'', '_scaled', '_LTSanalysis', '_burstanalysis'};
toCopyDir = {['/vtraces_old', num2str(versionNumber), '/'], ...
          ['/vtraces_scaled_old', num2str(versionNumber), '/'], ...
          ['/LTSanalysis_old', num2str(versionNumber), '/'], ...
          ['/burstanalysis_old', num2str(versionNumber), '/']};

%% Specify the prefix for subdirectories
prefix = 'EXAMPLE';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find path of matfile to use
% fullmatfilepath = fullfile(inFolder, filetouse);
fullmatfilepath = fullfile(inFolder, filetouse);
m = matfile(fullmatfilepath);

%% Extract info needed
fnrow = m.fnrow;
ltspeaktime = m.ltspeaktime;
spikesperpeak = m.spikesperpeak;
peakclass = m.peakclass;
peak2ndder = m.peak2ndder;

%% Find all traces with spikes but not LTSs (those spikes are presumably spontaneous spikes, but could also be noise)
All_Spontaneous_Spikes = intersect(find(spikesperpeak > 0), find(isnan(ltspeaktime)));
create_subdir_copy_files(All_Spontaneous_Spikes, fnrow, prefix, toCopySuffix, toCopyDir, inFolder, outFolder);

%% Find all traces with prominence lower than threshold (peakclass == 1) but with narrowness greater than threshold
Not_LTS_by_prom = intersect(find(peakclass == 1), find(peak2ndder < lts_thr));
create_subdir_copy_files(Not_LTS_by_prom, fnrow, prefix, toCopySuffix, toCopyDir, inFolder, outFolder);

%% Find all traces with LTS onset time > 3500 ms
Long_Latency_LTS = find(ltspeaktime > 3500);
create_subdir_copy_files(Long_Latency_LTS, fnrow, prefix, toCopySuffix, toCopyDir, inFolder, outFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% toCopyDir = {'/vtraces/', '/vtraces_scaled/', '/LTSanalysis/', '/burstanalysis/'};
% filetouse = 'dclampdatalog_take4.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%