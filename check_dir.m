function check_dir (directories)
%% Checks if needed directory(ies) exist and create them if not
% Usage: check_dir (directories)
% Arguments:
%       TODO
%
% Used by:    
%       /media/adamX/m3ha/data_dclamp/take4/find_initial_slopes.m
%
% File History:
% 2018-06-21 Modified from check_subdir.m
% TODO: Use print_or_show_message

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
% TODO

%% Check directory(ies)
if iscell(directories)
    for k = 1:numel(directories)
        if exist(directories{k}, 'dir') ~= 7
            mkdir(directories{k});
            fprintf('New subdirectory is made: %s\n\n', ...
                directories{k});
        end
    end
else
    if exist(directories, 'dir') ~= 7
        mkdir(directories);
        fprintf('New subdirectory is made: %s\n\n', ...
                directories);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
