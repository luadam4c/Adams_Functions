function indtofit = m3ha_find_ind_to_fit (fnrow, cellID, pharm, gincr, ...
                                            fitmode, infolder)
%% Find indices of fnrow in dclampdatalog_take4.mat that will be used for fitting
% Usage: indtofit = m3ha_find_ind_to_fit (fnrow, cellID, pharm, gincr, ...
%                                          fitmode, infolder)
% Arguments: 
%       fitmode     - 1 - all of g incr = 100%, 200%, 400%
%                   - 2 - all of g incr = 100%, 200%, 400% 
%                       but exclude cell-pharm-g_incr sets 
%                       containing problematic sweeps
%       infolder    - (opt) the directory that contains the special_cases folder
%                   must be a directory
%                   default == //media/adamX/m3ha/data_dclamp/take4/
%
% Used by:
%       cd/m3ha_dclampPassiveFitter.m
%       cd/m3ha_PlotHistogramsRefineThreshold.m
%       cd/m3ha_PlotCorrelations.m
%       cd/m3ha_dclampdatalog_analyze.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting4.m 
%           and later versions

% File History:
% 2016-09-06 Created
% 2016-09-09 Added _old
% 2016-09-13 Changed problematic sweeps
% 2016-09-13 Added fitmode
% 2016-10-15 Made infolder an optional argument; removed default for fitmode
% 2016-10-18 Skip infolder and related code if fitmode == 1
% 2016-10-19 Added Too_many_spontaneous_spikes and changed the name of folders
% 2016-10-20 Count each sweep only once (use only the _scaled file)
% 2017-01-16 Use all folders of the form TAKE_OUT_*/*.png
% 2017-05-22 Changed line width and indentation
% 2017-05-22 Renamed function FindIndToFit() -> find_ind_to_fit()
% 2018-10-01 infolder is now a relative path
% 2018-10-15 Renamed function find_ind_to_fit -> m3ha_find_ind_to_fit()
%               and moved to /home/Matlab/Adams_Functions
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 5
    error('A fitmode is required, type ''help m3ha_find_ind_to_fit'' for usage');
elseif isempty(fitmode) || ~isnumeric(fitmode) ...
        || ~(fitmode == 1 || fitmode == 2)
    error('fitmode out of range!, type ''help m3ha_find_ind_to_fit'' for usage');
elseif nargin >= 6 && ~isdir(infolder)
    error('infolder must be a directory!');
end

%% Set defaults for optional arguments
if fitmode == 2 && nargin < 6
    infolder = fullfile(pwd, 'take4/');
end

%% Extract information from inputs
ntraces = numel(fnrow);
if fitmode == 2
    infolderSp = fullfile(infolder, '/special_cases/');

    allFoldersToTakeOut = dir(fullfile(infolderSp, 'TAKE_OUT_*'));
    ct = 0;
    for folder = allFoldersToTakeOut'
        potentialFilesToTakeOut = ...
            dir(fullfile(infolderSp, folder.name, '*.png'));
        for file = potentialFilesToTakeOut'
            % Count each sweep only once (use only the _scaled file)
            if ~isempty(strfind(file.name, '_scaled')) ...
                && isempty(strfind(file.name, '_old'))        
                ct = ct + 1;
                matname = strrep(file.name, '_scaled.png', '.mat');
                filesToTakeOut{ct} = matname;
            end
        end
    end
end

%% Find indices corresponding to data for fitting
% Use only conductance amplitudes with 100%, 200% or 400% scaling 
%   (these are present in all experiments)
g_ind = find(gincr == 100 | gincr == 200 | gincr == 400);

% For fitmode == 2, find indices of sweeps that will not be used
if fitmode == 2
    nottouse_ind = [];
    for k = 1:ntraces
        if ~isempty(find(g_ind == k)) ...           % skip g incr of interest
            && isempty(find(nottouse_ind == k)) ... % don't recount same index
            && ismember(fnrow(k), filesToTakeOut)
            ntu_c_ind = find(cellID == cellID(k));
            ntu_p_ind = find(pharm == pharm(k));
            ntu_g_ind = find(gincr == gincr(k));
            ntu_cpg_ind = intersect(ntu_c_ind, intersect(ntu_p_ind, ntu_g_ind));    
            nottouse_ind = [nottouse_ind ntu_cpg_ind];
        end
    end
end

% Find indices to fit
if fitmode == 1
    indtofit = g_ind;
elseif fitmode == 2
    indtofit = setdiff(g_ind, nottouse_ind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%OLD CODE

%%% Looks_Like_Noise = {'D091710_0011_24.mat', 'E091710_0000_17.mat'};

% This condition was fixed 2016-09-12
Spontaneous_spikes_too_close_together_files = dir(fullfile(infolderSp, 'Spontaneous_spikes_too_close_together/*.png'));
ct = 0;
for file = Spontaneous_spikes_too_close_together_files'
    if isempty(strfind(file.name, '_LTSanalysis')) ...
        && isempty(strfind(file.name, '_burstanalysis')) ...
        && isempty(strfind(file.name, '_old'))
        ct = ct + 1;
        matname = strrep(file.name, '.png', '.mat');
        Spontaneous_spikes_too_close_together{ct} = matname;
    end
end

% This condition will be manually overridden 2016-09-13
% 20 sweeps remaining 2016-09-12
Spontaneous_spikes_that_did_not_fire_files = dir(fullfile(infolderSp, 'Spontaneous_spikes_that_did_not_fire/*.png'));
ct = 0;
for file = Spontaneous_spikes_that_did_not_fire_files'
    if isempty(strfind(file.name, '_LTSanalysis')) ...
        && isempty(strfind(file.name, '_burstanalysis')) ...
        && isempty(strfind(file.name, '_old'))
        ct = ct + 1;
        matname = strrep(file.name, '.png', '.mat');
        Spontaneous_spikes_that_did_not_fire{ct} = matname;
    end
end

% This condition will be manually overridden 2016-09-13
% 51 sweeps remaining 2016-09-13
Looks_like_noise_files = dir(fullfile(infolderSp, 'Looks_like_noise/*.png'));
ct = 0;
for file = Looks_like_noise_files'
    if isempty(strfind(file.name, '_LTSanalysis')) ...
        && isempty(strfind(file.name, '_burstanalysis')) ...
        && isempty(strfind(file.name, '_old'))
        ct = ct + 1;
        matname = strrep(file.name, '.png', '.mat');
        Looks_like_noise{ct} = matname;
    end
end

    % 2 sweeps remaining 2016-10-19
    Spontaneous_LTSs_or_bursts_fixed_files = dir(fullfile(infolderSp, 'Spontaneous_LTSs_or_bursts_fixed/*.png'));
    ct = 0;
    for file = Spontaneous_LTSs_or_bursts_fixed_files'
        if isempty(strfind(file.name, '_LTSanalysis')) ...
            && isempty(strfind(file.name, '_burstanalysis')) ...
            && isempty(strfind(file.name, '_old'))        % ignore duplicate files that were made for comparison
            ct = ct + 1;
            matname = strrep(file.name, '.png', '.mat');
            Spontaneous_LTSs_or_bursts_fixed{ct} = matname;
        end
    end

    % Problematic sweeps for which the entire set of cell-pharm-g_incr condition will be left out
    % 23 sweeps remaining 2016-09-13
    Noisy_recording_files = dir(fullfile(infolderSp, 'CONTESTED_TAKE_OUT_Noisy_recording/*.png'));
    ct = 0;
    for file = Noisy_recording_files'
        if ~isempty(strfind(file.name, '_scaled')) ...
            && isempty(strfind(file.name, '_old'))        % Count each sweep only once (use only the _scaled file)
            ct = ct + 1;
            matname = strrep(file.name, '_scaled.png', '.mat');
            Noisy_recording{ct} = matname;
        end
    end

    % 54 sweeps remaining 2016-10-19
    Spontaneous_LTSs_or_bursts_files = dir(fullfile(infolderSp, 'CONTESTED_TAKE_OUT_Spontaneous_LTSs_or_bursts/*.png'));
    ct = 0;
    for file = Spontaneous_LTSs_or_bursts_files'
        if ~isempty(strfind(file.name, '_scaled')) ...
            && isempty(strfind(file.name, '_old'))        % Count each sweep only once (use only the _scaled file)
            ct = ct + 1;
            matname = strrep(file.name, '_scaled.png', '.mat');
            Spontaneous_LTSs_or_bursts{ct} = matname;
        end
    end

    % 122 sweeps remaining 2016-10-19
    Too_many_spontaneous_spikes_files = dir(fullfile(infolderSp, 'CONTESTED_TAKE_OUT_Too_many_spontaneous_spikes/*.png'));
    ct = 0;
    for file = Too_many_spontaneous_spikes_files'
        if ~isempty(strfind(file.name, '_scaled')) ...
            && isempty(strfind(file.name, '_old'))        % Count each sweep only once (use only the _scaled file)
            ct = ct + 1;
            matname = strrep(file.name, '_scaled.png', '.mat');
            Too_many_spontaneous_spikes{ct} = matname;
        end
    end

 ...
            || ismember(fnrow(k), Spontaneous_LTSs_or_bursts_files) ...
            || ismember(fnrow(k), Too_many_spontaneous_spikes_files))

    infolder = '//media/adamX/m3ha/data_dclamp/take4/';

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
