function check_subdir (outfolder, subdir)
%% Checks if needed subdirectory(ies) exist in outfolder
% Usage: check_subdir (outfolder, subdir)
% Arguments:
%		%%% TO DO
%
% Used by:	
%		/media/adamX/m3ha/data_dclamp/take4/find_special_cases.m
%		/media/adamX/m3ha/optimizer4gabab/optimizer_4compgabab.m
%		/home/Matlab/Adams_Functions/find_passive_params.m
%		/home/Matlab/Adams_Functions/find_istart.m
%		/home/Matlab/Adams_Functions/find_IPSC_peak.m
%		/home/Matlab/Adams_Functions/find_LTS.m
%		/home/Matlab/Adams_Functions/create_subdir_copy_files.m
%       /home/Matlab/minEASE.m
%
% File History:
% 2016-11-02 Created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
%%% TO DO

%% Check subdirectory(ies)
if iscell(subdir)
	for k = 1:numel(subdir)
		if exist(fullfile(outfolder, subdir{k}), 'dir') ~= 7
			mkdir(fullfile(outfolder, subdir{k}));
			fprintf('New subdirectory is made: %s\n\n', ...
				fullfile(outfolder, subdir{k}));
		end
	end
else
	if exist(fullfile(outfolder, subdir), 'dir') ~= 7
		mkdir(fullfile(outfolder, subdir));
		fprintf('New subdirectory is made: %s\n\n', ...
				fullfile(outfolder, subdir));
	end
end

