function create_subdir_copy_files (setIndices, fileNames, prefix, toCopySuffix, toCopyDir, inFolder, outFolder)
%% Creates subdirectories and copies figure files
% Usage: create_subdir_copy_files (setIndices, fileNames, prefix, toCopySuffix, toCopyDir, inFolder, outFolder)
% Arguments:
%       TODO
%
% Requires:    
%       cd/check_subdir.m
%
% Used by:    
%       cd/m3ha_find_special_cases.m

% File History:
% 2016-12-06 Moved from /media/adamX/m3ha/data_dclamp/find_special_cases.m
% TODO: Input parser
% TODO: Make more general by adding FileType (default 'png') 
%           as an optional parameter-value pair

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO
setName = inputname(1);

%% TODO
subDirName = strcat(prefix, '_', setName);

check_subdir(outFolder, subDirName);

%% TODO
for i = 1:length(setIndices)
    for k = 1:numel(toCopySuffix)
        pngFileNames = strrep(fileNames{setIndices(i)}, ...
                    '.mat', [toCopySuffix{k}, '.png']);
        pngPathsFrom = fullfile(inFolder, toCopyDir{k}, pngFileNames);
        pngPathsTo = fullfile(outFolder, subDirName, pngFileNames);
        if ~isfile(pngPathsTo)
            copyfile(pngPathsFrom, pngPathsTo);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%