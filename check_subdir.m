function check_subdir (parentDirectory, subDirectories)
%% Checks if needed subdirectory(ies) exist in parentDirectory
% Usage: check_subdir (parentDirectory, subDirectories)
% Arguments:
%      TODO
%
% Used by:    
%       /media/adamX/m3ha/data_dclamp/take4/find_special_cases.m
%       /media/adamX/m3ha/optimizer4gabab/optimizer_4compgabab.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting22.m
%       /home/Matlab/Adams_Functions/find_passive_params.m
%       /home/Matlab/Adams_Functions/find_istart.m
%       /home/Matlab/Adams_Functions/find_IPSC_peak.m
%       /home/Matlab/Adams_Functions/find_LTS.m
%       /home/Matlab/Adams_Functions/create_subdir_copy_files.m
%       /home/Matlab/Adams_Functions/create_input_file.m
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/EEG_gui/EEG_gui.m
%       /home/Matlab/function_template.m
%
% File History:
% 2016-11-02 Created
% 2018-06-19 Changed tabs to spaces
% TODO: Use print_or_show_message

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
% TODO

%% Check subdirectory(ies)
if iscell(subDirectories)
    for k = 1:numel(subDirectories)
        if exist(fullfile(parentDirectory, subDirectories{k}), 'dir') ~= 7
            mkdir(fullfile(parentDirectory, subDirectories{k}));
            fprintf('New subdirectory is made: %s\n\n', ...
                fullfile(parentDirectory, subDirectories{k}));
        end
    end
else
    if exist(fullfile(parentDirectory, subDirectories), 'dir') ~= 7
        mkdir(fullfile(parentDirectory, subDirectories));
        fprintf('New subdirectory is made: %s\n\n', ...
                fullfile(parentDirectory, subDirectories));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
