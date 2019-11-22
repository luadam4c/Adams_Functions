function m3ha_phase_plane_analysis(datafilename, xlimits, isLTS)
%% Plots phase planes from an m3ha .mat file
% Arguments: 
%       'D091710_0007_15.mat', [1000 8910]
%
% Used by:
%       cd/m3ha_test_phase_plane_analysis.m

% File History:
% 2016-07-30 Created

mafw = 30;	% width in ms for the moving average filter

data_dir = '//media/adamX/m3ha/data_dclamp/take4/matfiles/';
outfolder = '//media/adamX/m3ha/data_dclamp/take4/phase_plots/';
if nargin > 2
	if isLTS	
		outfolder = '//media/adamX/m3ha/data_dclamp/take4/phase_plots/true_LTS';
	else	
		outfolder = '//media/adamX/m3ha/data_dclamp/take4/phase_plots/false_LTS';
	end
end

if nargin < 1
	error('No filename specified');
else
	data_num = strrep(datafilename, '.mat', '');
	m = matfile(fullfile(data_dir, datafilename));
end

sims = m.d_orig(2, 1) - m.d_orig(1, 1);
if nargin < 2
	nd_mfrs = size(m.d_mfrs, 1);
	xlimits = [m.d_mfrs(2, 1) m.d_mfrs(nd_mfrs, 1)];
end
ind = round(xlimits(1)/sims):round(xlimits(2)/sims);

ndp = length(ind);
ndp_mafw = round(mafw/sims);
ndp2 = ndp - 2 * ndp_mafw;

tvec0 = m.d_orig(ind, 1);
vvec0 = m.d_orig(ind, 4);
vvec1 = m.d_mf(ind, 4);

dvvec0 = diff(vvec0)./diff(tvec0);
ddvvec0 = diff(dvvec0)./diff(tvec0(2:end));

dvvec1 = diff(vvec1)./diff(tvec0);
dvvec1_sm = smooth(dvvec1, ndp_mafw);
dvvec1_sm_cut = dvvec1_sm(2 + ndp_mafw:ndp - ndp_mafw);
ddvvec1 = diff(dvvec1_sm_cut)./diff(tvec0(2 + ndp_mafw:ndp - ndp_mafw));

vvec3 = smooth(vvec1, ndp_mafw);
dvvec3 = diff(vvec3)./diff(tvec0);
dvvec3_sm = smooth(dvvec3, ndp_mafw);
dvvec3_sm_cut = dvvec3_sm(2 + ndp_mafw:ndp - ndp_mafw);
ddvvec3 = diff(dvvec3_sm_cut)./diff(tvec0(2 + ndp_mafw:ndp - ndp_mafw));
np2der = min(ddvvec3);
fprintf('np2der = %10.5f\n', np2der);

h = figure; 
set(h, 'Visible', 'off');
set(h, 'Name', 'Phase plane analysis');
clf(h);
subplot(3,1,1)
plot(vvec0(2:end), dvvec0, 'b-', 'LineWidth', 0.5);
ylim([-2 2]);
xlabel('Voltage (mV)');
ylabel('Voltage slope (V/s)');
title('Phase plane analysis, original trace')
subplot(3,1,2)
plot(vvec1(2:end), dvvec1, 'g-', 'LineWidth', 0.5);
ylim([-2 2]);
xlabel('Voltage (mV)');
ylabel('Voltage slope (V/s)');
title('Phase plane analysis, median-filtered trace')
subplot(3,1,3)
plot(vvec3(2:end), dvvec3, 'r-', 'LineWidth', 0.5);
ylim([-1 1]);
xlabel('Voltage (mV)');
ylabel('Voltage slope (V/s)');
title('Phase plane analysis, moving-average-filtered trace')
figname = fullfile(outfolder, [data_num, '_phase', '.png']);
saveas(h, figname);

h = figure;
set(h, 'Visible', 'off');
set(h, 'Name', 'LTS peak analysis, original trace');
clf(h);
subplot(3,1,1)
plot(tvec0, vvec0, 'b-', 'LineWidth', 0.5); hold on
xlim(xlimits);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('LTS peak analysis, original trace')
subplot(3,1,2)
plot(tvec0(2:end), dvvec0, 'k-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('dV/dT');
subplot(3,1,3)
plot(tvec0(3:end), ddvvec0, 'k-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('d2V/dT2');
figname = fullfile(outfolder, [data_num, '_orig', '.png']);
saveas(h, figname);

h = figure;
set(h, 'Visible', 'off');
set(h, 'Name', 'LTS peak analysis, median-filtered trace');
clf(h);
subplot(3,1,1)
plot(tvec0, vvec1, 'g-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('LTS peak analysis, median-filtered trace')
subplot(3,1,2)
plot(tvec0(2:end), dvvec1, 'k-', 'LineWidth', 0.5);
plot(tvec0(2:end), dvvec1_sm, 'r-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('dV/dT');
subplot(3,1,3)
plot(tvec0(3 + ndp_mafw:ndp - ndp_mafw), ddvvec1, 'k-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('d2V/dT2');
figname = fullfile(outfolder, [data_num, '_mf', '.png']);
saveas(h, figname);

h = figure;
set(h, 'Visible', 'off');
set(h, 'Name', 'LTS peak analysis, moving-average-filtered trace');
clf(h);
subplot(3,1,1)
plot(tvec0, vvec3, 'r-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('LTS peak analysis, moving-average-filtered trace')
subplot(3,1,2)
plot(tvec0(2:end), dvvec3, 'k-', 'LineWidth', 0.5);
plot(tvec0(2:end), dvvec3_sm, 'r-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('dV/dT');
subplot(3,1,3)
plot(tvec0(3 + ndp_mafw:ndp - ndp_mafw), ddvvec3, 'k-', 'LineWidth', 0.5);
xlim(xlimits);
xlabel('Time (ms)');
ylabel('d2V/dT2');
figname = fullfile(outfolder, [data_num, '_maf', '.png']);
saveas(h, figname);

display(['This file is finished: ', datafilename]);

