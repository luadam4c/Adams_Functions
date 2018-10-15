function [IPSC_offset_old, IPSC_offset_old2, IPSC_ind_old, IPSC_ind_old2] = find_istart_old (tvec0, gvec1s, ivec0s, istartwin, plotflag, vvec0s, outfolder, filebase)
%% Finds time of current application from a series of current vectors
% Usage: [IPSC_offset_old, IPSC_offset_old2, IPSC_ind_old, IPSC_ind_old2] = find_istart_old (tvec0, gvec1s, ivec0s, istartwin, plotflag, vvec0s, outfolder, filebase)
% Arguments:	
%
% Used by:
%		/media/adamX/m3ha/data_dclamp/dclampDataExtractor.m
%
% 2016-11-07 Moved from dclampDataExtractor.m

mafw1 = 5;	% width in ms for the moving average filter for finding IPSC offsets
slth = 0.1;	% slope threshold for finding IPSC offset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sims = tvec0(2) - tvec0(1);
nswps = size(ivec0s, 2);

% Find and plot IPSC offsets
ndp_mafw1 = round(mafw1/sims);
gstdvec = std(gvec1s, 0, 2);			% Compute the standard deviation over all sweeps
gstdvec_sm = smooth(gstdvec, ndp_mafw1);	% Smooth with a moving average filter spanning 5 ms
gstdvec_sm_slope = diff(gstdvec_sm)/sims;	% Slope of smoothed standard deviation
istdvec = std(ivec0s, 0, 2);			% Compute the standard deviation over all sweeps
istdvec_sm = smooth(istdvec, ndp_mafw1);	% Smooth with a moving average filter spanning 5 ms
istdvec_sm_slope = diff(istdvec_sm)/sims;	% Slope of smoothed standard deviation
ind = round(istartwin(1)/sims):round(istartwin(2)/sims);	% indices of interest
g_first_rise_pt = find(gstdvec_sm_slope(ind) > slth, 1);	% This is the old method for comparison
if isempty(g_first_rise_pt)
	IPSC_offset_old = NaN;
	IPSC_ind_old = round(1000/sims);
else
	IPSC_offset_old = (g_first_rise_pt - 1)*sims; 	% IPSC offset in ms, using 0.1 as the slope threshold
	IPSC_ind_old = round((istartwin(1) + IPSC_offset_old)/sims);	% index of IPSC application
end
i_first_dip_pt = find(istdvec_sm_slope(ind) > slth, 1);		% Slope of standard deviation is positive
if isempty(i_first_dip_pt)
	IPSC_offset_old2 = NaN;				% To detect in .csv file
	IPSC_ind_old2 = round(1000/sims);			% For the rest of the analysis to work smoothly
else
	IPSC_offset_old2 = (i_first_dip_pt - 1)*sims; 	% IPSC offset in ms, using 0.1 as the slope threshold
	IPSC_ind_old2 = round((istartwin(1) + IPSC_offset_old2)/sims);	% index of IPSC application
end

if plotflag
	h = figure(3000);
	set(h, 'Visible', 'off');
	set(h, 'Name', 'IPSC offset analysis (obsolete)');
	clf(h);
	subplot(3,1,1);
	for swp = 1:nswps					% Plot each median-filtered conductance trace
		plot(tvec0(ind), gvec1s(ind, swp)); hold on; 
	end
	plot(tvec0(ind), gstdvec(ind), 'LineWidth', 2, 'Color', 'r');
	plot(tvec0(ind), gstdvec_sm(ind), 'LineWidth', 2, 'Color', 'g');
	plot(tvec0(ind), gstdvec_sm_slope(ind), 'LineWidth', 2, 'Color', 'b');
	if ~isempty(i_first_dip_pt)
		plot(tvec0(IPSC_ind_old2), gstdvec_sm(IPSC_ind_old2), ...
		    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
		plot(tvec0(IPSC_ind_old), gstdvec_sm(IPSC_ind_old), ...
		    'LineWidth', 2, 'Color', 'k', 'Marker', '+', 'MarkerSize', 12);
	end
	title(['IPSC offset analysis (obsolete) for ', filebase]);
	xlabel('Time (ms)')
	ylabel('Conductance (uS)')
	xlim(istartwin);
	ylim([-5 10]);
	subplot(3,1,2);
	for swp = 1:nswps					% Plot each current trace
		plot(tvec0(ind), ivec0s(ind, swp)); hold on; 
	end
	plot(tvec0(ind), istdvec(ind), 'LineWidth', 2, 'Color', 'r');
	plot(tvec0(ind), istdvec_sm(ind), 'LineWidth', 2, 'Color', 'g');
	plot(tvec0(ind), istdvec_sm_slope(ind), 'LineWidth', 2, 'Color', 'b');
	if ~isempty(i_first_dip_pt)
		plot(tvec0(IPSC_ind_old2), istdvec_sm(IPSC_ind_old2), ...
		    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
		plot(tvec0(IPSC_ind_old), istdvec_sm(IPSC_ind_old), ...
		    'LineWidth', 2, 'Color', 'k', 'Marker', '+', 'MarkerSize', 12);
	end
	xlabel('Time (ms)')
	ylabel('Current (pA)')
	xlim(istartwin);
	ylim([-100 10]);				
	subplot(3,1,3);
	for swp = 1:nswps					% Plot each voltage trace
		plot(tvec0(ind), vvec0s(ind, swp)); hold on; 
	end
	if ~isempty(i_first_dip_pt)
		plot(tvec0(IPSC_ind_old2), vvec0s(IPSC_ind_old2, 1), ...
		    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
		plot(tvec0(IPSC_ind_old), vvec0s(IPSC_ind_old, 1), ...
		    'LineWidth', 2, 'Color', 'k', 'Marker', '+', 'MarkerSize', 12);
	end
	xlabel('Time (ms)')
	ylabel('Voltage (mV)')
	xlim(istartwin);
	ylim([-90 -50]);				
	figname = fullfile(outfolder, '/IPSCoffset_old/', [filebase, '_IPSCoffsetold', '.png']);
	saveas(h, figname);
	close(h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%