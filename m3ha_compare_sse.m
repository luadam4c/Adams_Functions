% m3ha_compare_sse.m
%% Compares the sum-of-squares error between the different trace correction strategies with boxplots after the outliers are taken out (with remove_outliers.m)
%
% Requires:	
%		cd/remove_outliers.m

% File History:
% 2016-12-08 Created
% 2018-06-11 Updated usage of remove_outliers.m
%


%% Set parameters
wl2IQR = 10;				% the ratio of whisker length to interquartile range

%% Set label for each condition
condition_labels = {'Rs, Gdata', 'Rs, Gtheo', 'Voff, Gdata', 'Voff, Gtheo'};

%% Set folders for reading and saving files
infolder = '//media/adamX/m3ha/data_dclamp/take4/backup/';
outfolder = '//media/adamX/m3ha/data_dclamp/take4/boxplots_all/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(outfolder, 'dir') ~= 7
	mkdir(outfolder);
end

%% Extract cell ID #s
m3 = matfile(fullfile(infolder, 'dclampdatalog_take4_old13.mat'));
cellidrow = m3.cellidrow;

%% Find vectors to compare
m1 = matfile(fullfile(infolder, 'trace_comparison_old13-1.mat'));
m2 = matfile(fullfile(infolder, 'trace_comparison_old13-5.mat'));
val_all = [m1.Rs1', m1.Rs2', m2.Voff1', m2.Voff2'];
sse_all = [m1.sse1', m1.sse2', m2.sse1', m2.sse2'];

%% Remove outliers of SSE (possibly something wrong with trace)
[sse_all_trimmed, indices_left] = ...
	remove_outliers (sse_all, 'OutlierMethod', 'boxplot', 'WL2IQR', wl2IQR);
cellidrow = cellidrow(indices_left);
val_all_trimmed = val_all(indices_left, :);
means = repmat(mean(val_all_trimmed, 1), size(val_all_trimmed, 1), 1);
zscore_all_trimmed = (val_all_trimmed - means) ./ means;

%% Make box plot for total sse across conditions
h = figure(floor((1+rand())*10^6));
clf(h);
boxplot(sse_all_trimmed, 'Labels', condition_labels);
ylabel('Sum-of-squares error');
suptitle(['SSE with outliers greater than ', num2str(wl2IQR), 'x the IQR taken out']);
figname = fullfile(outfolder, 'compare_sse_all.png');
saveas(h, figname);

%% Make box plots for sse across cells
parfor c = 1:numel(condition_labels)
	h = figure(floor((1+rand())*10^6));
	clf(h);
	boxplot(sse_all_trimmed(:, c)', cellidrow, 'Plotstyle', 'compact');
	ylabel('Sum-of-squares error');
	suptitle(['SSE with outliers greater than ', num2str(wl2IQR), 'x the IQR taken out; ', condition_labels{c}]);
	figname = fullfile(outfolder, ['compare_sse_', strrep(condition_labels{c}, ', ', '_') '_c.png']);
	saveas(h, figname);
end

%% Make box plot for Rs or Voff across conditions
h = figure(floor((1+rand())*10^6));
clf(h);
boxplot(zscore_all_trimmed, 'Labels', condition_labels);
ylabel('z score');
suptitle(['Rs or Voff with outliers greater than ', num2str(wl2IQR), 'x the IQR taken out']);
figname = fullfile(outfolder, 'compare_Rs_Voff_all.png');
saveas(h, figname);

%% Make box plots for Rs or Voff across cells
parfor c = 1:numel(condition_labels)
	h = figure(floor((1+rand())*10^6));
	clf(h);
	boxplot(val_all_trimmed(:, c)', cellidrow, 'Plotstyle', 'compact');
	if c == 1 || c == 2
		ylabel('Series Resistance (MOhm)');
		suptitle([condition_labels{c}, ' with outliers greater than ', num2str(wl2IQR), 'x the IQR taken out']);
	elseif c == 3 || c == 4
		ylabel('Voltage offset (mV)');
		suptitle([condition_labels{c}, ' with outliers greater than ', num2str(wl2IQR), 'x the IQR taken out']);
	end
	figname = fullfile(outfolder, ['compare_', strrep(condition_labels{c}, ', ', '_') '_c.png']);
	saveas(h, figname);
end

%{
% OLD CODE

h = findobj(gca, 'tag', 'Outliers');
delete(h);

[val_all_trimmed, indices_left1] = remove_outliers (val_all, wl2IQR);
cellidrow1 = cellidrow(indices_left1);
cellidrow2 = cellidrow(indices_left2);
	boxplot(sse_all_trimmed(:, c)', cellidrow2, 'Plotstyle', 'compact');
	boxplot(val_all_trimmed(:, c)', cellidrow1, 'Plotstyle', 'compact');

%}
