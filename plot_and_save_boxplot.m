function plot_and_save_boxplot (h, vec, veclabel, outfolder, filename, title_mod, group_vec, groupn_label)
%% Plots a box plot from a grouped vector according to group
% Usage: plot_and_save_boxplot (h, vec, veclabel, outfolder, filename, title_mod, group_vec, groupn_label)
%
% Requires:	
%
% Used by:	
%		/media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
% 
% 2016-12-08 Adapted from plot_and_save_histogram.m
% 

%% Check arguments
%%% UNFINISHED

%% Set defaults
%%% UNFINISHED
if nargin < 3
	veclabel = 'something';
end
if nargin < 6
	title_mod = '';
end
if nargin < 7
	group_vec = [];
end
if nargin < 8
	group_labels = [];
end

%% Create figure and plot boxplot
% histg.m is from /home/Matlab/Adams_Functions/
h = figure(h);
set(h, 'Name', ['Distribution of ', veclabel]);
clf(h);
boxplot(vec, group_vec, 'Plotstyle', 'compact');			% plots boxplot compact style
%legend('Location', 'eastoutside');
text(1, -0.1, groupn_label, 'Units', 'normalized');			% group number label
ylabel(veclabel);							% data vector label
title(['Comparison of ', veclabel, title_mod]);

if nargin >= 4
	figname = fullfile(outfolder, filename);
	saveas(h, figname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%