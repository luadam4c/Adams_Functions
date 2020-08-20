function CreateSubdirAndCopyPNGFiles (setindices, fnrow, prefix, tocopy_suffix, tocopy_dir, infolder, outfolder)
%% Create subdirectory and copy figure files
% Usage: CreateSubdirAndCopyPNGFiles (setindices, fnrow, prefix, tocopy_suffix, tocopy_dir, infolder, outfolder)
% Arguments: 	%%% TO DO
%
% Requires:	
%		cd/check_subdir.m
%
% Used by:	
%		/media/adamX/m3ha/data_dclamp/find_special_cases.m
%
% 2016-12-06 Moved from /media/adamX/m3ha/data_dclamp/find_special_cases.m
% TODO: Make more general by passing FileType ('png') as optional argument


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setname = inputname(1);
for i = 1:length(setindices)
	for k = 1:numel(tocopy_suffix)
		pngfile = strrep(fnrow{setindices(i)}, '.mat', [tocopy_suffix{k}, '.png']);
		pngfilepath = fullfile(infolder, tocopy_dir{k}, pngfile);
		subdir = sprintf('/%s_%s/', prefix, setname);
		check_subdir(outfolder, subdir);
		pngfilepath_cp = fullfile(outfolder, subdir, pngfile);
		if exist(pngfilepath_cp) ~= 2
			copyfile(pngfilepath, pngfilepath_cp);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

