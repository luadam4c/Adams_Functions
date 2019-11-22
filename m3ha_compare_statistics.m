function m3ha_compare_statistics(new_vn, old_vn)
% Used to find and copy traces with new stats such as peak classifications, LTS peak times and spikes per peak
% Usage: m3ha_compare_statistics(new_vn, old_vn)
%%% TODO: Please test m3ha_compare_statistics(14, 13)
% Arguments:	new_vn	- new version number
%			must be an integer between 1 and max_vn
%		old_vn	- old version number
%			must be an integer between 1 and new_vn
%
% Requires:	
%		cd/find_in_strings.m

% File History:
% 2016-09-11 Created
% 2016-10-14 Added version number
% 2016-10-17 Removed the need to adjust for ioffset_old(k) (for version number 4 and above)
% 2016-10-17 Added the need to adjust for the change in istart_ind definition (for version number 5 and above)
% 2016-12-05 Added peakclass as a comparison stat
% 2016-12-05 Made the comparison dependent on filename instead of index
% 2017-03-23 BT - Made new_vn, old_vn arguments, read newer versions from //media/adamX/m3ha/data_dclamp/take4/backup/, compares new peakclass against equivalent old peakclass
% 2017-04-04 AL - Fixed the names of the outfolders
% 

max_vn = 14;			% maximum version number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
	error('No new version number specified!');
elseif ~isnumeric(new_vn) || mod(new_vn, 1) ~= 0
	error('New version number must be an integer!');
end;
if nargin < 2
	error('No old version number specified!');
elseif ~isnumeric(old_vn) || mod(old_vn, 1) ~= 0
	error('Old version number must be an integer!');
end 
if new_vn < 1 || new_vn > max_vn
	error('New version number must be between 1 and max version number!');
end
if old_vn < 1 || old_vn > new_vn
	error('Old version number must be between 1 and new version number!');
end

%% Set folders for reading and saving files
infolder1 = '//media/adamX/m3ha/data_dclamp/take4/';
infolder2 = '//media/adamX/m3ha/data_dclamp/take4/backup/'; 
infolder3 = ['//media/adamX/m3ha/data_dclamp/take4/backup/vtraces_old', num2str(new_vn), '/'];
infolder4 = ['//media/adamX/m3ha/data_dclamp/take4/backup/vtraces_old', num2str(old_vn), '/'];
infolder5 = ['//media/adamX/m3ha/data_dclamp/take4/backup/LTSanalysis_old', num2str(new_vn), '/'];
infolder6 = ['//media/adamX/m3ha/data_dclamp/take4/backup/LTSanalysis_old', num2str(old_vn), '/'];
infolder7 = ['//media/adamX/m3ha/data_dclamp/take4/backup/burstanalysis_old', num2str(new_vn), '/'];
infolder8 = ['//media/adamX/m3ha/data_dclamp/take4/backup/burstanalysis_old', num2str(old_vn), '/'];
infolder9 = ['//media/adamX/m3ha/data_dclamp/take4/backup/vtraces_scaled_old', num2str(new_vn), '/'];
infolder10 = ['//media/adamX/m3ha/data_dclamp/take4/backup/vtraces_scaled_old', num2str(old_vn), '/'];
outfolder1 = ['//media/adamX/m3ha/data_dclamp/take4/to_compare_lts_old', ...
			num2str(new_vn), 'against_old', num2str(old_vn), '/'];
outfolder2 = ['//media/adamX/m3ha/data_dclamp/take4/to_compare_spp_old', ...
			num2str(new_vn), 'against_old', num2str(old_vn), '/'];
outfolder3 = ['//media/adamX/m3ha/data_dclamp/take4/to_compare_pkclass_old', ...
			num2str(new_vn), '_against_old', num2str(old_vn), '/'];
if exist(outfolder1, 'dir') ~= 7
	mkdir(outfolder1);
end
if exist(outfolder2, 'dir') ~= 7
	mkdir(outfolder2);
end
if exist(outfolder3, 'dir') ~= 7
	mkdir(outfolder3);
end

%% Set files for comparing
new_file = fullfile(infolder2, ['dclampdatalog_take4_old', num2str(new_vn), '.mat']);
old_file = fullfile(infolder2, ['dclampdatalog_take4_old', num2str(old_vn), '.mat']);

%% Load data
d1 = matfile(new_file);
d2 = matfile(old_file);
ntraces = numel(d1.fnrow);

%% Load filenames
fnrow_new = d1.fnrow;
fnrow_old = d2.fnrow;
ltspeaktime_new = d1.ltspeaktime;
ltspeaktime_old = d2.ltspeaktime;
spikesperpeak_new = d1.spikesperpeak;
spikesperpeak_old = d2.spikesperpeak;
peakclass_new = d1.peakclass;
peakclass_old = d2.peakclass;

%% Version difference to change in peakclass?
version_diff = old_vn <= 13 && new_vn >= 14 && (range(peakclass_new) - range(peakclass_old)) >= 2;

%% Find differences in LTS peak times
ct = 0;
for k = 1:ntraces
	filename = fnrow_new{k};
	ltst_new = ltspeaktime_new(k);
	ltst_old = ltspeaktime_old(find_in_strings(filename, fnrow_old));
	% Only use if comparing against versions 1~3
	%	if (ltst_new ~= ltst_old ...
	%		&& abs(ltst_new - ltst_old - d2.ioffset_old(k)) > 0.1) ...
	% Only use if comparing against version 4~5
	if (~isnan(ltst_new) && ~isnan(ltst_old) && (ltst_new - ltst_old) > 0.099) ...
		|| (isnan(ltst_new) && ~isnan(ltst_old)) ...
		|| (~isnan(ltst_new) && isnan(ltst_old)) 
	% To use in the future
	% if (~isnan(ltst_new) && ~isnan(ltst_old) && ltst_new ~= ltst_old) ...
		ct = ct + 1;
		diff_lts{ct} = filename;		% record file name
		diff_lts_val(ct, 1) = ltst_old;		% record old ltst value
		diff_lts_val(ct, 2) = ltst_new;		% record new ltst value
		fprintf('%s\n', filename);
		pngfile = strrep(filename, '.mat', '.png');
		pngfile_new = fullfile(infolder3, pngfile);
		pngfile_old = fullfile(infolder4, pngfile);
		pngfile_new_cp = fullfile(outfolder1, pngfile);
		pngfile_old_cp = fullfile(outfolder1, strrep(filename, '.mat', '_old.png'));
		copyfile(pngfile_new, pngfile_new_cp);
		copyfile(pngfile_old, pngfile_old_cp);
	end
end
fprintf('There are %d files with LTS peak times different from before!\n', ct);

%% Find differences in spikes per peak
ct = 0;
for k = 1:ntraces
	filename = fnrow_new{k};
	spp_new = spikesperpeak_new(k);
	spp_old = spikesperpeak_old(find_in_strings(filename, fnrow_old));
	if spp_new ~= spp_old	% spikes per peak is always a number
		ct = ct + 1;
		diff_spp{ct} = filename;
		diff_spp_val(ct, 1) = spp_old;
		diff_spp_val(ct, 2) = spp_new;
		fprintf('%s\n', filename);
		pngfile = strrep(filename, '.mat', '.png');
		pngfile_new = fullfile(infolder3, pngfile);
		pngfile_old = fullfile(infolder4, pngfile);
		pngfile_new_cp = fullfile(outfolder2, pngfile);
		pngfile_old_cp = fullfile(outfolder2, strrep(filename, '.mat', '_old.png'));
		copyfile(pngfile_new, pngfile_new_cp);
		copyfile(pngfile_old, pngfile_old_cp);
		pngfile = strrep(filename, '.mat', '_LTSanalysis.png');
		pngfile_new = fullfile(infolder5, pngfile);
		pngfile_old = fullfile(infolder6, pngfile);
		pngfile_new_cp = fullfile(outfolder2, pngfile);
		pngfile_old_cp = fullfile(outfolder2, strrep(filename, '.mat', '_LTSanalysis_old.png'));
		copyfile(pngfile_new, pngfile_new_cp);
		copyfile(pngfile_old, pngfile_old_cp);
		pngfile = strrep(filename, '.mat', '_burstanalysis.png');
		pngfile_new = fullfile(infolder7, pngfile);
		pngfile_old = fullfile(infolder8, pngfile);
		pngfile_new_cp = fullfile(outfolder2, pngfile);
		pngfile_old_cp = fullfile(outfolder2, strrep(filename, '.mat', '_burstanalysis_old.png'));
		copyfile(pngfile_new, pngfile_new_cp);
		copyfile(pngfile_old, pngfile_old_cp);
		pngfile = strrep(filename, '.mat', '_scaled.png');
		pngfile_new = fullfile(infolder9, pngfile);
		pngfile_old = fullfile(infolder10, pngfile);
		pngfile_new_cp = fullfile(outfolder2, pngfile);
		pngfile_old_cp = fullfile(outfolder2, strrep(filename, '.mat', '_scaled_old.png'));
		copyfile(pngfile_new, pngfile_new_cp);
		copyfile(pngfile_old, pngfile_old_cp);	
	end
end
fprintf('There are %d files with spikes per peak different from before!\n', ct);

%% Find differences in peakclass
ct = 0;
for k = 1:ntraces
	filename = fnrow_new{k};
	pkclass_new = peakclass_new(k);
	pkclass_old = peakclass_old(find_in_strings(filename, fnrow_old));
	if version_diff && pkclass_old >= 4 && pkclass_old <= 7
		pkclass_old = pkclass_old + 2;
	end
	if pkclass_new ~= pkclass_old	% spikes per peak is always a number
		ct = ct + 1;
		diff_pkclass{ct} = filename;
		diff_pkclass_val(ct, 1) = pkclass_old;
		diff_pkclass_val(ct, 2) = pkclass_new;
		fprintf('%s\n', filename);
		pngfile = strrep(filename, '.mat', '.png');
		pngfile_new = fullfile(infolder3, pngfile);
		pngfile_old = fullfile(infolder4, pngfile);
		pngfile_new_cp = fullfile(outfolder3, pngfile);
		pngfile_old_cp = fullfile(outfolder3, strrep(filename, '.mat', '_old.png'));
		copyfile(pngfile_new, pngfile_new_cp);
		copyfile(pngfile_old, pngfile_old_cp);
	end
end
fprintf('There are %d files with peakclass different from before!\n', ct);


%{
OLD CODE:

function m3ha_compare_statistics(varargin)
% Usage: m3ha_compare_statistics(varargin)

%% Set version number to compare against
vn = 13;			% version number

infolder1 = '//media/adamX/m3ha/data_dclamp/take4/';
infolder2 = '//media/adamX/m3ha/data_dclamp/take4/backup/'; 
infolder3 = '//media/adamX/m3ha/data_dclamp/take4/vtraces/';
infolder4 = ['//media/adamX/m3ha/data_dclamp/take4/backup/vtraces_old', num2str(vn), '/'];
infolder5 = '//media/adamX/m3ha/data_dclamp/take4/LTSanalysis/';
infolder6 = ['//media/adamX/m3ha/data_dclamp/take4/backup/LTSanalysis_old', num2str(vn), '/'];
infolder7 = '//media/adamX/m3ha/data_dclamp/take4/burstanalysis/';
infolder8 = ['//media/adamX/m3ha/data_dclamp/take4/backup/burstanalysis_old', num2str(vn), '/'];
infolder9 = '//media/adamX/m3ha/data_dclamp/take4/vtraces_scaled/';
infolder10 = ['//media/adamX/m3ha/data_dclamp/take4/backup/vtraces_scaled_old', num2str(vn), '/'];
outfolder1 = ['//media/adamX/m3ha/data_dclamp/take4/to_compare_lts_new_against_old', num2str(vn), '/'];
outfolder2 = ['//media/adamX/m3ha/data_dclamp/take4/to_compare_spp_new_against_old', num2str(vn), '/'];
outfolder3 = ['//media/adamX/m3ha/data_dclamp/take4/to_compare_pkclass_new_against_old', num2str(vn), '/'];

new_file = fullfile(infolder1, 'dclampdatalog_take4.mat');
%}
