function [ndps, sims, tvec0, gvec0s, ivec0s, vvec0s, gvec1s, ivec1s, vvec1s, ndps2, sims2, tvec2, gvec2s, ivec2s, vvec2s, vvec3s] = load_matfiles_part (newinfolder, datfn, nswps)
% Load set of matfiles and return relevant info
% Usage: [ndps, sims, tvec0, gvec0s, ivec0s, vvec0s, gvec1s, ivec1s, vvec1s, ndps2, sims2, tvec2, gvec2s, ivec2s, vvec2s, vvec3s] = load_matfiles_part (newinfolder, datfn, nswps)
%
% Used by:
%		/media/adamX/m3ha/data_dclamp/dclampDataExtractor.m
%
% 2016-11-07 Moved from dclampDataExtractor.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check whether matfile exists
thisfile = fullfile(newinfolder, [datfn, '_1.mat']);
if exist(thisfile, 'file') ~= 2
	error(['This mat file: ', thisfile, ' is missing!!']);
end

% Load matfile
m = matfile(thisfile);

% Obtain time vectors from .mat data of first sweep
tvec0 = m.d_orig(:, 1);
tvec2 = m.d_mfrs(:, 1);
ndps = length(tvec0);
ndps2 = length(tvec2);
sims = tvec0(2) - tvec0(1);
sims2 = tvec2(2) - tvec2(1);

% Load relevant .mat data for all sweeps

gvec0s = zeros(ndps, nswps);
ivec0s = zeros(ndps, nswps);
vvec0s = zeros(ndps, nswps);
gvec1s = zeros(ndps, nswps);
ivec1s = zeros(ndps, nswps);
vvec1s = zeros(ndps, nswps);
gvec2s = zeros(ndps2, nswps);
ivec2s = zeros(ndps2, nswps);
vvec2s = zeros(ndps2, nswps);
vvec3s = zeros(ndps, nswps);			
parfor swp = 1:nswps		% FOR each sweep
	filebase = [datfn, '_', num2str(swp)];
	thisfile = fullfile(newinfolder, [filebase, '.mat']);
	if ~exist(thisfile, 'file')
		error(['This mat file: ', thisfile, ' is missing!!']);
	end
	m = matfile(thisfile);

	% Extract data
	gvec0s(:, swp) = m.d_orig(:, 2);
	ivec0s(:, swp) = m.d_orig(:, 3);
	vvec0s(:, swp) = m.d_orig(:, 4);
	gvec1s(:, swp) = m.d_mf(:, 2);
	ivec1s(:, swp) = m.d_mf(:, 3);
	vvec1s(:, swp) = m.d_mf(:, 4);
	gvec2s(:, swp) = m.d_mfrs(:, 2);
	ivec2s(:, swp) = m.d_mfrs(:, 3);
	vvec2s(:, swp) = m.d_mfrs(:, 4);
	vvec3s(:, swp) = m.d_mfmaf(:, 4);
end

