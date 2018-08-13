function [tvec0, gvec0, ivec0, vvec0, tvec1, gvec1, ivec1, vvec1, tvec2, gvec2, ivec2, vvec2, tvec3, gvec3, ivec3, vvec3] = plot_traces_mat (filename, xlimits)
%% Plot traces from mat file
% Usage: [tvec0, gvec0, ivec0, vvec0, tvec1, gvec1, ivec1, vvec1, tvec2, gvec2, ivec2, vvec2, tvec3, gvec3, ivec3, vvec3] = plot_traces_mat (filename, xlimits)
% Arguments:	filename	- either full path or within current directory, 
%					e.g. '/media/adamX/m3ha/data_dclamp/take4/matfiles/D091710_0007_15.mat'
%		xlimits		- (opt) time interval to show, e.g., [1800 2000]
%
% Used by:
%		/media/adamX/m3ha/data_dclamp/test_sweep.m

% 20160729 - Created
% 20160909 - Added mfmaf
% 20160922 - removed data_dir, modified filename requirements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check for required arguments
if nargin < 1
	error('No filename specified!');
elseif ~ischar(filename)
	error('Filename must be a char array in single quotes!');
end

%% Create full path to abf file robustly
[filepath, filebase, ~] = fileparts(filename);
matfilename = strcat(filebase, '.mat');
if isempty(filepath)
	filepath = pwd;
end
filename_full = fullfile(filepath, matfilename);
if exist(filename_full) ~= 2
	error('This mat file doesn''t exist!');
end
fprintf('Full path to mat file: %s\n', filename_full);

%% Load data
m = matfile(filename_full);

tvec0 = m.d_orig(:, 1);
gvec0 = m.d_orig(:, 2);
ivec0 = m.d_orig(:, 3);
vvec0 = m.d_orig(:, 4);
tvec1 = m.d_mf(:, 1);
gvec1 = m.d_mf(:, 2);
ivec1 = m.d_mf(:, 3);
vvec1 = m.d_mf(:, 4);
tvec2 = m.d_mfrs(:, 1);
gvec2 = m.d_mfrs(:, 2);
ivec2 = m.d_mfrs(:, 3);
vvec2 = m.d_mfrs(:, 4);
tvec3 = m.d_mfmaf(:, 1);
gvec3 = m.d_mfmaf(:, 2);
ivec3 = m.d_mfmaf(:, 3);
vvec3 = m.d_mfmaf(:, 4);

sims = tvec0(2) - tvec0(1);
sims2 = tvec2(2) - tvec2(1);
if nargin < 2
	xlimits = [tvec2(1) tvec2(end)];
end
ind = round(xlimits(1)/sims):round(xlimits(2)/sims);
ind2 = round(xlimits(1)/sims2):round(xlimits(2)/sims2);

h = figure(1);
set(h, 'Name', 'Conductance traces');
clf(h);
plot(tvec0(ind), gvec0(ind), 'b-', 'LineWidth', 0.5); hold on
%plot(tvec0(ind), gvec1(ind), 'g-', 'LineWidth', 0.5);
plot(tvec2(ind2), gvec2(ind2), 'g.', 'LineWidth', 0.5);
plot(tvec0(ind), gvec3(ind), 'r-', 'LineWidth', 0.5);
title('Conductance traces');
%legend('raw trace', 'median filtered (width 2.5 ms)', 'median filtered (width 2.5 ms) and resampled', 'median-filtered then moving-average-filtered', 'Location', 'SouthOutside')
legend('raw trace', 'median filtered (width 2.5 ms) and resampled', 'median-filtered then moving-average-filtered', 'Location', 'SouthOutside')
xlim(xlimits);
xlabel('Time (ms)');
ylabel('Conductance (uS)');

h = figure(2);
set(h, 'Name', 'Current traces');
clf(h);
plot(tvec0(ind), ivec0(ind), 'b-', 'LineWidth', 0.5); hold on
%plot(tvec0(ind), ivec1(ind), 'g.-', 'LineWidth', 0.5);
plot(tvec2(ind2), ivec2(ind2), 'g.-', 'LineWidth', 0.5);
plot(tvec0(ind), ivec3(ind), 'r-', 'LineWidth', 0.5);
title('Current traces');
%legend('raw trace', 'median filtered (width 10 ms)', 'median filtered (width 10 ms) and resampled', 'median-filtered then moving-average-filtered', 'Location', 'SouthOutside')
legend('raw trace', 'median filtered (width 10 ms) and resampled', 'median-filtered then moving-average-filtered', 'Location', 'SouthOutside')
xlim(xlimits);
xlabel('Time (ms)');
ylabel('Current (nA)');

h = figure(3);
set(h, 'Name', 'Voltage traces');
clf(h);
plot(tvec0(ind), vvec0(ind), 'b-', 'LineWidth', 0.5); hold on
%plot(tvec0(ind), vvec1(ind), 'g-', 'LineWidth', 0.5);
plot(tvec2(ind2), vvec2(ind2), 'g.', 'LineWidth', 0.5);
plot(tvec0(ind), vvec3(ind), 'r-', 'LineWidth', 0.5);
title('Voltage traces');
%legend('raw trace', 'median filtered (width 30 ms)', 'median filtered (width 30 ms) and resampled', 'median-filtered then moving-average-filtered', 'Location', 'SouthOutside')
legend('raw trace', 'median filtered (width 30 ms) and resampled', 'median-filtered then moving-average-filtered', 'Location', 'SouthOutside')
xlim(xlimits);
xlabel('Time (ms)');
ylabel('Voltage (mV)');

h = figure(4);
vvec4 = -vvec3;
vvec5 = flipud(vvec3);
vvec6 = -flipud(vvec3);
set(h, 'Name', 'LTS analysis');
clf(h);
subplot(2, 2, 1)
plot(tvec0(ind), vvec3(ind), 'r-', 'LineWidth', 0.5);
xlim([4000 4500]);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('original');
subplot(2, 2, 2)
plot(tvec0(ind), vvec4(ind), 'r-', 'LineWidth', 0.5);
xlim([4000 4500]);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('negative');
subplot(2, 2, 3)
plot(tvec0(ind), vvec5(ind), 'r-', 'LineWidth', 0.5);
xlim([tvec0(end)-4500 tvec0(end)-4000]);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('flipped');
subplot(2, 2, 4)
plot(tvec0(ind), vvec6(ind), 'r-', 'LineWidth', 0.5);
xlim([tvec0(end)-4500 tvec0(end)-4000]);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('negative of flipped');

