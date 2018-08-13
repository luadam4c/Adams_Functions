function plot_and_save_histogram (h, vec, veclabel, countlabel, outfolder, filename, nbins, title_mod, class_vec, class_labels)
%% Plots and saves a stacked histogram for a vector and color code according to class
% Usage: plot_and_save_histogram (h, vec, veclabel, countlabel, outfolder, filename, nbins, title_mod, class_vec, class_labels)
%
% Requires:	
%		/home/Matlab/Adams_Functions/histg.m
% Used by:	
%		/media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
% 
% 2016-??-?? Created
% 2016-10-14 Moved to /home/Matlab/Adams_Functions/
% 2016-12-01 If classes doesn't exist, plot regular histogram instead of stacked histogram
% 2016-12-01 Added countlabel
% 2016-12-04 Removed legend for regular histograms
% 2016-12-08 Changed argument from fitmode to title_mod
% 

%% Check arguments
% UNFINISHED

%% Set defaults
% UNFINISHED
if nargin < 3
	veclabel = 'something';
end
if nargin < 7
	nbins = 50;
end
if nargin < 8
	title_mod = '';
end
if nargin < 9
	class_vec = [];
end
if nargin < 10
	class_labels = [];
end

%% Create figure and plot histogram
% histg.m is from /home/Matlab/Adams_Functions/
h = figure(h);
set(h, 'Name', ['Distribution of ', veclabel]);
clf(h);
if ~isempty(class_vec) && ~isempty(class_labels)
	opt.nbins = nbins;			% needed for histg
	opt.group_names = class_labels';	% needed for histg
	histg(vec, class_vec, opt);		% plots stacked histogram
	legend('Location', 'eastoutside');
else
	histogram(vec, nbins);			% plots regular histogram
end
xlabel(veclabel)
ylabel(countlabel)
title(['Distribution of ', veclabel, title_mod]);

if nargin >= 6
	figname = fullfile(outfolder, filename);
	saveas(h, figname);
end

% Old codes
%{

if fitmode == 0
	title(['Distribution of ', veclabel, ' (all)']);
elseif fitmode == 1
	title(['Distribution of ', veclabel, ' (100%, 200%, 400% g incr)']);
elseif fitmode == 2
	title(['Distribution of ', veclabel, ' (for fitting)']);
end

%}
