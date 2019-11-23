function [cpa_ap, g_sc, i_sc] = m3ha_resave_sweeps (abffullfn, nswps_supposed, ljp, IPSC_start_time, gabab_amp, gabab_Trise, gabab_TfallFast, gabab_TfallSlow, gabab_w, cpmid, outfolder)
%% Extracts .abf data and resave as .mat file in the m3ha format
% Usage: [cpa_ap, g_sc, i_sc] = m3ha_resave_sweeps (abffullfn, nswps_supposed, ljp, IPSC_start_time, gabab_amp, gabab_Trise, gabab_TfallFast, gabab_TfallSlow, gabab_w, cpmid, outfolder)
% Arguments:
%	abffullfn	- full file path to abf file
%	nswps_cpv	- 
% 	ljp		- liquid junction potential used to correct voltage trace (mV)
%	IPSC_start_time - 
%	gabab_amp	- 
%	gabab_Trise	- 
%	gabab_TfallFast - 
%	gabab_TfallSlow - 
%	gabab_w		- 
%	cpmid		- 
%	outfolder	- 
%
% Requires:
%		~/Downloaded_Functions/abf2load.m
%		cd/compute_gabab_conductance.m
%		cd/rescale_vec.m
%
% Used by:
%		cd/m3ha_parse_dclamp_data.m
%
% 2016-11-07 Moved from m3ha_parse_dclamp_data.m
% 2016-12-13 Reversed sign of LJP

%% Parameters used for data reorganization
rsims = 1;	% resampling interval in ms (1 kHz)
mfw1 = 2.5;	% width in ms for the median filter for PClamp noise (conductance traces)
mfw2 = 10;	% width in ms for the median filter for corrupted data (current traces)
mfw3 = 30;	% width in ms for the median filter for spikes (voltage traces)
mafw2 = 30;	% width in ms for the moving average filter for finding narrowest voltage peaks
mvw = 0.5;	% width in ms for calculating mean voltage
precision = 0.2;	% precision of the scaling factor for correcting conductance traces

%% Problematic abf files
flipped_files = {'F092210_0000', 'F092210_0001', 'F092210_0002', ...
			'F092210_0003', 'F092210_0004', 'F092210_0005', ...
			'F092210_0006', 'F092210_0007'};	% The current and conductance channels were flipped

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract name of set and number of sweeps used
[~, datfn, ~] = fileparts(abffullfn);

% Load .abf data
[d, si, ~] = abf2load(abffullfn);	% CLK - these were pclamp10 files.
			% CLK - sacrifice use of "ep" epoch info extracter in abfload2.m
ndps = size(d, 1);				% Total number of data points in each sweep
nswps = size(d, 3);				% Total number of sweeps in this abf file
if nswps ~= nswps_supposed
	error('Number of sweeps is incorrect!!');
end

% Reorganize data
sims = si/1000;					% sampling interval in milliseconds (0.1 ms)
tvec0 = sims*(1:size(d, 1))';			% time vector from 0.1 ms to 9000.0 or 9500.0 ms
tvec2 = rsims*(1:round(tvec0(end)/rsims))';	% time vector from 1 ms to 9000 or 9500 ms
ndps2 = length(tvec2);
gvec0 = zeros(ndps, nswps);
ivec0 = zeros(ndps, nswps);
vvec0 = zeros(ndps, nswps);
gvec1 = zeros(ndps, nswps);
ivec1 = zeros(ndps, nswps);
vvec1 = zeros(ndps, nswps);
gvec2 = zeros(ndps2, nswps);
ivec2 = zeros(ndps2, nswps);
vvec2 = zeros(ndps2, nswps);
gvec3 = zeros(ndps, nswps);
ivec3 = zeros(ndps, nswps);
vvec3 = zeros(ndps, nswps);
cpa_ap = zeros(1, nswps);			% approximate current pulse amplitude (pA)
g_sc = ones(1, nswps);				% scaling factor for conductance traces, default is 1
i_sc = ones(1, nswps);				% scaling factor for current traces, default is 1
parfor swp = 1:nswps				% FOR each sweep
	% Record current and voltage traces; 
	% NOTE: some data points seem corrupted
	if ismember(datfn, flipped_files)
		gvec0_temp = d(:, 2, swp);	% conductance trace in pS OR nS
		ivec0_temp = d(:, 3, swp);	% current trace in pA
	else
		gvec0_temp = d(:, 3, swp);	% conductance trace in pS OR nS
		ivec0_temp = d(:, 2, swp);	% current trace in pA
	end

	% Convert conductance traces to nS if necessary
	if max(gvec0_temp) > 100
		gvec0_temp = gvec0_temp / 1000;	% convert from pS to nS
		g_sc(swp) = 0.001;
	end

	% Fix conductance traces that were arbitrarily scaled by comparing with theoretical trace
	gvec0_theo = compute_gabab_conductance(tvec0, IPSC_start_time, ...
			gabab_amp(swp), gabab_Trise(swp), ...
			gabab_TfallFast(swp), gabab_TfallSlow(swp), gabab_w(swp));
	gvec0_temp(gvec0_temp < 0) = 0;		% remove negative values from gvec0_temp first
	[gvec0(:, swp), new_sc] = rescale_vec(gvec0_temp, gvec0_theo, precision);
	g_sc(swp) = g_sc(swp) * new_sc;
	
	% Fix current traces that were arbitrarily scaled
	cpmidind = round(cpmid/sims):round((cpmid + mvw)/sims);		% indices for measuring cp peak
	cpa_ap(swp) = mean(ivec0_temp(cpmidind));	% approximate cpa
	if cpa_ap(swp) < -105 && cpa_ap(swp) > -145	% about -125 pA
		i_sc(swp) = 0.4;			% was arbitrarily scaled by 2.5
	elseif cpa_ap(swp) < -230 && cpa_ap(swp) > -270  % about -250 pA
		i_sc(swp) = 0.2;			% was arbitrarily scaled by 5
	end
	ivec0(:, swp) = ivec0_temp * i_sc(swp);

	% Do LJP correction (-10 mV) for voltage traces
	vvec0(:, swp) = d(:, 1, swp) - ljp;	% voltage trace LJP-corrected

	% Median filter conductance traces to get rid of PClamp noise
	gvec1(:, swp) = medfilt1(gvec0(:, swp), round(mfw1/sims));

	% Median filter current traces to get rid of corrupted data
	ivec1(:, swp) = medfilt1(ivec0(:, swp), round(mfw2/sims));

	% Median filter voltage traces to get rid of spikes
	vvec1(:, swp) = medfilt1(vvec0(:, swp), round(mfw3/sims));

	% Resample all traces at 1 kHz for fitting use (Christine's traces)
	gvec2(:, swp) = interp1(tvec0, gvec1(:, swp), tvec2, 'linear');
	ivec2(:, swp) = interp1(tvec0, ivec1(:, swp), tvec2, 'linear');
	vvec2(:, swp) = interp1(tvec0, vvec1(:, swp), tvec2, 'linear');

	% Moving-average-filter median-filtered traces for taking derivatives
	gvec3(:, swp) = smooth(gvec1(:, swp), round(mafw2/sims));
	ivec3(:, swp) = smooth(ivec1(:, swp), round(mafw2/sims));
	vvec3(:, swp) = smooth(vvec1(:, swp), round(mafw2/sims));
end

% Save data as .mat file
for swp = 1:nswps		% FOR each sweep
	filebase = [datfn, '_', num2str(swp)];
	d_orig = zeros(ndps, 4);
	d_mf = zeros(ndps, 4);
	d_mfrs = zeros(ndps2, 4);
	d_mfmaf = zeros(ndps, 4);
	d_orig(:, 1) = tvec0;
	d_orig(:, 2) = gvec0(:, swp);
	d_orig(:, 3) = ivec0(:, swp);
	d_orig(:, 4) = vvec0(:, swp);
	d_mf(:, 1) = tvec0;
	d_mf(:, 2) = gvec1(:, swp);
	d_mf(:, 3) = ivec1(:, swp);
	d_mf(:, 4) = vvec1(:, swp);
	d_mfrs(:, 1) = tvec2;
	d_mfrs(:, 2) = gvec2(:, swp);
	d_mfrs(:, 3) = ivec2(:, swp);
	d_mfrs(:, 4) = vvec2(:, swp);
	d_mfmaf(:, 1) = tvec0;
	d_mfmaf(:, 2) = gvec3(:, swp);
	d_mfmaf(:, 3) = ivec3(:, swp);
	d_mfmaf(:, 4) = vvec3(:, swp);
	filename = [filebase, '.mat'];
	newdatafn = fullfile(outfolder, '/matfiles/', filename);
	save(newdatafn, 'd_orig', 'd_mf', 'd_mfrs', 'd_mfmaf', ...
			'cpa_ap', 'g_sc', 'i_sc', '-v7.3');
				% Cannot use parfor if save is in the loop
end

%{


%}
