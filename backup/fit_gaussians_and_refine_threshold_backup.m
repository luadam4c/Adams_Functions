function [BestModel, numComponents, Mu_best, Stdev_best, Prop_best, minAIC, new_thr, alt_thr, new_peakclass, min_threshold] = fit_gaussians_and_refine_threshold (data, max_numComponents, old_thr, thr_min, label, outfolder, filename, fitmode, peakclass, peakclass_labels, min_threshold)
%% Fits data to Gaussian mixture models and finds the optimal number of components
% arguments: data must be a column vector of values; 
% 		max_numComponents is defaulted to be 5
%		the rest are only necessary if the histogram is to be plotted and saved
%		Examine cell plot must have figname rmse_*R/F*_row_Fit_traces and threshold plot must have figname rmse_*R/F*_row_Fit_threshold
% 
% Requires:	
%		/home/Matlab/Adams_Functions/histg.m
%		/home/Matlab/Brians_Functions/reorder.m
% Used by:	
%		/media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%
% 2016-08-03 Created
% 2016-09-15 Added fitmode
% 2016-09-30 Added peakclass & peakclass_labels
% 2016-09-30 Took out spontaneous spikes 
% 2016-11-01 Replaced ProbDistUnivParam with makedist, also from the Statistics and Machine Learning Toolbox
% 2017-02-08 - BT - Changed legend location for RMSE graph fitting
% 2017-02-22 - BT - Identified RMSEs above threshold, returns modified peakclass for RMSE
% 2017-04-29 Don't apply reorder.m unless plotting RMSE
% 2017-05-01 - BT - remove_cells for cells to be examined between RMSE plots. 
% 2017-05-05 - BT - Examine cell plot must have figname rmse_*R/F*_row_Fit_traces and threshold plot must have figname rmse_*R/F*_row_Fit_threshold
% TODO: add input parser and make arguments except data parameter-value pairs with sensible default values


%% Set parameters
prec = 10^-4;				% Precision
remove_cells = [5 14 19];		% Cell #s to examine in RMSE plotting

%% Set up log file
logfile = fullfile(outfolder, strrep(filename, '.png', ['_fit_gaussians_', num2str(max_numComponents), 'maxcomplog.txt']));
fid = fopen(logfile, 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
	functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
	functionsdirectory = '/scratch/al4ng/Matlab/';
else
	error('Valid functionsdirectory does not exist!');
end
addpath(fullfile(functionsdirectory, '/Brians_Functions/'));		% for reorder.m

%% Take out spontaneous spikes (peak class == 3)
ind_nospont = find(peakclass ~= 3);
data = data(ind_nospont);
peakclass = peakclass(ind_nospont);
ind_inf = find(isinf(data));
data = data(~isinf(data));
ind_nospont(ind_inf) = [];
peakclass(ind_inf) = [];

%% Extract data info
npts = length(data);

%% Find largest order of magnitude in data
% Set dmax and dmin to 99th and 1st percentile to match histg
dmax = prctile(data, 99);
dmin = prctile(data, 1);
dabs_max = max([abs(dmax) abs(dmin)]);
ordomag = floor(log10(dabs_max));

%% Find proper bounds and bins for histogram
binw = 10^(ordomag - 1);		% Bin width
left = floor(dmin/binw)*binw;		% left bound
right = ceil(dmax/binw)*binw;		% right bound
edges = (left:binw:right)';

%% Find area of histogram
harea = binw * npts;

%% Find the best Gaussian mixture fit
if nargin < 2
	max_numComponents = 5;
end
AIC = zeros(max_numComponents,1);
GMModels = cell(max_numComponents,1);			% Preallocation of Gaussian mixture models
Mus = cell(max_numComponents,1);			% Preallocation of means of Gaussian mixture components
Stdevs = cell(max_numComponents,1);			% Preallocation of standard deviations of Gaussian mixture components
Props = cell(max_numComponents,1);			% Preallocation of proportions of Gaussian mixture components
options = statset('MaxIter', 1000);
rng(1);					% For reproducibility
for k = 1:max_numComponents
	GMModels{k} = fitgmdist(data, k, 'Options', options);	
	Mu = zeros(k, 1);
	Stdev = zeros(k, 1);
	Prop = zeros(k, 1);
	for i = 1:k
		Mu(i) = GMModels{k}.mu(i);
		Stdev(i) = sqrt(GMModels{k}.Sigma(1,1,i));
	Prop(i) = GMModels{k}.ComponentProportion(1,i);
	end
	[Mus{k}, I_s] = sort(Mu);
	Stdevs{k} = Stdev(I_s);
	Props{k} = Prop(I_s);
	fprintf(fid, '\n GM Mean(s), Standard deviation(s) & Proportions(s) for %i Component(s)\n', k);
	for i = 1:k
		fprintf(fid, '%g\t%g\t%g\n', Mus{k}(i), Stdevs{k}(i), Props{k}(i));
	end
	AIC(k)= GMModels{k}.AIC;
end
[minAIC, numComponents] = min(AIC);
BestModel = GMModels{numComponents};
Mu_best = Mus{numComponents};
Stdev_best = Stdevs{numComponents};
Prop_best = Props{numComponents};
fprintf(fid, '\n GM Model selected has %d components \n', numComponents);
fprintf(fid, '\n Mean(s)\tStandard deviation(s)\tProportions(s)\n');
for i = 1:numComponents
	fprintf(fid, '%g\t%g\t%g\n', Mu_best(i), Stdev_best(i), Prop_best(i));
end

%% Extract info about the Gaussian mixture fit
x = (left:prec:right)';
Pdf = pdf(BestModel, x);
Normal = cell(numComponents, 1);
Normalpdf = cell(numComponents, 1);
for i = 1:numComponents
	Normal{i} = makedist('Normal', 'mu', Mu_best(i), 'sigma', Stdev_best(i));	% create a normal distribution for this component
	Normalpdf{i} = pdf(Normal{i}, x);						% create the probability density function
end

%% Find new thresholds
if numComponents == 3					% 2 possible thresholds

	% Find the minimum between lower two Gaussians
	left1 = max(floor(Mu_best(1)/prec)*prec, thr_min);	% left bound; threshold can't be lower than thr_min
	right1 = max(ceil(Mu_best(2)/prec)*prec, left1);	% right bound
	xx1 = (left1:prec:right1)';
	Pdf_part1 = pdf(BestModel, xx1);
	[minprob1 minprob1_ind] = min(Pdf_part1);
	new_thr1 = xx1(minprob1_ind);

	% Find the minimum between upper two Gaussians
	left2 = max(floor(Mu_best(2)/prec)*prec, thr_min);	% left bound; threshold can't be lower than thr_min
	right2 = max(ceil(Mu_best(3)/prec)*prec, left2);	% right bound
	xx2 = (left2:prec:right2)';
	Pdf_part2 = pdf(BestModel, xx2);
	[minprob2 minprob2_ind] = min(Pdf_part2);
	new_thr2 = xx2(minprob2_ind);

	% Use the closest value to the old threshold as the new threshold
	[~, choice] = min(abs([new_thr1 new_thr2] - old_thr));
	if choice == 1
		new_thr = new_thr1;
		alt_thr = new_thr2;
	elseif choice == 2
		new_thr = new_thr2;
		alt_thr = new_thr1;
	else
		error('threshold is wrong!')
	end

	fprintf(fid, 'New LTS threshold is %g V^2/s^2\n', new_thr);
	fprintf(fid, 'Alternate LTS threshold is %g V^2/s^2\n', alt_thr);
	fclose(fid);
 
elseif numComponents == 4				% 3 possible thresholds
else
	comp1 = find(Mu_best < old_thr, 1, 'last');	% Gaussian peak left of old_thr
	comp2 = find(Mu_best > old_thr, 1);		% Gaussian peak right of old_thr
	value1 = Mu_best(comp1);
	value2 = Mu_best(comp2);

	left1 = max(floor(value1/prec)*prec, thr_min);	% left bound; threshold can't be lower than thr_min
	right1 = ceil(value2/prec)*prec;		% right bound
	xx = (left1:prec:right1)';
	Pdf_part = pdf(BestModel, xx);
	[minprob minprob_ind] = min(Pdf_part);
	new_thr = xx(minprob_ind);
	fprintf(fid, 'New LTS threshold is %g V^2/s^2\n', new_thr);
	fclose(fid);
end

%% Plot and save histogram with fitted Gaussian mixture pdf
if fitmode == 0
	fitsuffix = ' (all)';
elseif fitmode == 1
	fitsuffix = ' (100%, 200%, 400% g incr)';
elseif fitmode == 2
	fitsuffix = ' (for fitting)';
end

new_peakclass = peakclass;	% copy peakclass to new_peakclass as template
new_peakclass_labels = peakclass_labels;	% copy peakclass_labels to new_peakclass_labels as template if not changed later

if nargin > 3		%%% TODO: Why 3? BT - in the original source
	h = figure(400);
	set(h, 'Visible', 'on');
	set(h, 'Name', ['Distribution of ', label]);
	clf(h)
	opt.edges = edges;			% needed for histg
	opt.group_names = peakclass_labels';	% needed for histg
	if length(peakclass_labels) >= 10 && strcmp(label, 'RMSE (mV) in the falling phase');	% RMSE falling phase plot analysis plots
		if strfind(filename, 'rmse_F_row_Fit_threshold')		% marker for plotting above/below threshold 
			new_peakclass_labels = {'Above Threshold', 'Below Threshold'};	% labels for above/below threshold plot
			opt.group_names = new_peakclass_labels';	% for histg.m
		elseif strfind(filename, 'rmse_F_row_Fit_traces') && size(remove_cells,2) > 0 	% plotting specific cells
			new_peakclass_labels = {'All other cells'};		% Make labels for cells to remove
			for x = 1:size(remove_cells,2)
				new_peakclass_labels{x+1} = ['Cell #' num2str(remove_cells(x))];
			end
			opt.group_names = new_peakclass_labels';	% for histg.m
		end
	end
	histg(data, new_peakclass, opt);		% plots stacked histogram
	for i = 1:numComponents				% plots gaussian components
		plot(x, Normalpdf{i} * Prop_best(i) * harea, 'w', 'Displayname', ['component #', num2str(i)]);
	end
	full_dist = plot(x, Pdf * harea, 'r', 'LineWidth', 1, 'Displayname', 'full distribution');	% plots full distribution
	if numComponents == 3
        	line([new_thr new_thr], [0 npts], 'Color', 'r', 'LineStyle', '--', ...
            		'Displayname', ['New threshold = ', num2str(new_thr)]); hold on;
		line([alt_thr alt_thr], [0 npts], 'Color', 'b', 'LineStyle', '--', ...
			'Displayname', ['Alternate threshold = ', num2str(alt_thr)]);
	end
	new_thr = -1; alt_thr = -1;              % temp		%%% TODO: annotate here
	xlim([left right]); 
	ylim([0 max(Pdf * harea)]); 
	xlabel(label)
	ylabel('# of sweeps')
	title(['Distribution of ', label, fitsuffix]);

	if length(peakclass_labels) < 10		% Legend for 2nd derivative, using peakclass
		legend('Location', 'northwest'); 
  	elseif length(peakclass_labels) >= 10 && ...
		strcmp(label, 'RMSE (mV) in the rising phase')
					% Legend for rising RMSE, location northwest overlaps plot; peakclass = cellidrow
		[min_val, min_ind] = min(full_dist.YData);	% Minimum y value of data to determine threshold
		x_min = full_dist.XData(min_ind);	% Corresponding x value of minimum y value
		if strfind(filename, 'rmse_R_row_Fit_threshold')	% above/below threshold RMSE plot
			new_peakclass(data >= x_min) = 2;	% Change peakclasses to visibly identify traces above threshold
			new_peakclass(data < x_min) = 1;	% All data above and below min x value altered
			new_peakclass_labels = {'Above Threshold', 'Below Threshold'};
		elseif strfind(filename, 'rmse_R_row_Fit_traces') && size(remove_cells,2) > 0	% cells are being examined between rising/falling plots
			new_peakclass_labels = {'All other cells'};		% Make labels for cells to remove
			not_remove_cells_indices = {};				% indices of cells not in remove_cells
			for x = 1:size(remove_cells,2)
				new_peakclass_labels{x+1} = ['Cell #' num2str(remove_cells(x))];	% make cell labels
				not_remove_cells_indices{x} = find(peakclass ~= remove_cells(x));	
								% track indices of cells not being examined
			end
			opt.group_names = new_peakclass_labels';	% for histg.m
			new_peakclass(intersectm(not_remove_cells_indices)) = 1;	
								% change all indices of cells not being examined to 1
			for x = 1:size(remove_cells,2)		% cells being examined given new classes for histg.m
				new_peakclass(peakclass == remove_cells(x)) = x+1;
			end
		end
		clf(h);	% Clear figure to replot with new peakclasses
		histg(data, new_peakclass, opt);	% Replot with new peakclass identifier
        	new_peakclass = new_peakclass';		% Inverted to pass as argument from rising to falling phase
		for i = 1:numComponents			% Plot each Gaussian component
			plot(x, Normalpdf{i} * Prop_best(i) * harea, 'w', 'Displayname', ['component #', num2str(i)]);
        	end
        	full_dist = plot(x, Pdf * harea, 'r', 'LineWidth', 1, 'Displayname', 'full distribution');	
							% Plot the final fit
		xlim([left right]); 
		ylim([0 max(Pdf * harea)]); 
		xlabel(label)
		ylabel('# of sweeps')
		title(['Distribution of ', label, fitsuffix]);
		legend(new_peakclass_labels, 'Location', 'eastoutside');	% Revised legend, original if no analysis criteria met
		line([x_min x_min], [0 npts], 'Color', 'r', 'LineStyle', '--', ...
		    'Displayname', ['New threshold = ', num2str(new_thr)]); hold on;	% Plot threshold line at minimum in rising
		text(1.1, 0.95, ['Rising Threshold: ' num2str(x_min)], 'Units', 'normalized');		
										% Text at what RMSE the threshold is at
    		text(1.1, 0.85, ['Components: ' num2str(max_numComponents)], 'Units', 'normalized');	
										% Text how many components are present
		min_threshold = x_min;	% Pass as argument to falling phase for labelling
  	elseif length(peakclass_labels) >= 10 && ...
		strcmp(label, 'RMSE (mV) in the falling phase')		% Legend for falling RMSE without replotting
		legend(new_peakclass_labels, 'Location', 'eastoutside');	% Revised legend, original if no analysis criteria met
		text(1.1, 0.95, ['Rising Threshold: ' num2str(min_threshold)], 'Units', 'normalized');		% Text at what RMSE the threshold is at
		text(1.1, 0.85, ['Components: ' num2str(max_numComponents)], 'Units', 'normalized');	% Text how many components are present	
	end
end
if nargin > 4		%%% TODO: Why 4? BT - also in original source
	figname = fullfile(outfolder, filename);
	saveas(h, figname);
	xlim([-0.009 0]);
	ylim([0 500]);
	title(['Distribution of ', label, ', zoomed in', fitsuffix]);
	figname = fullfile(outfolder, strrep(filename, '.png', '_zoomed.png'));
	saveas(h, figname);
	close(h);
else
	close(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
%% Old Code

dmax = max(data); 
dmin = min(data);

	histogram(data, edges);

%% Group data by peak class
ind_cl = cell(1, 7);					% Indices for each peak class
data_cl = cell(1, 7);					% data for each peak class
color_cl = {'y', 'm', 'c', 'r', 'g', 'b', 'k'};
for pkcl = 1:7
	ind_cl{pkcl} = find(peakclass == pkcl);
	data_cl{pkcl} = data(ind_cl{pkcl});
end

	for pkcl = 1:7
		histogram(data_cl{pkcl}, edges, 'FaceColor', color_cl{pkcl}, 'EdgeColor', color_cl{pkcl});
	end

	if numComponents == 3
		legend(['New threshold = ', num2str(new_thr)], ...
			['Alternate threshold = ', num2str(alt_thr)], ...
			'Peak class #1', 'Peak class #2', 'Peak class #3', ...
			'Peak class #4', 'Peak class #5', 'Peak class #6', ...
			'Peak class #7', 'Location', 'northwest');
	else
		legend(['New threshold = ', num2str(new_thr)], ...
			'Peak class #1', 'Peak class #2', 'Peak class #3', ...
			'Peak class #4', 'Peak class #5', 'Peak class #6', ...
			'Peak class #7', 'Location', 'northwest');
	end

Normal{i} = ProbDistUnivParam('normal', [Mu_best(i) Stdev_best(i)]);
new_peakclass(peakclass ~= 5 & peakclass ~= 14 & peakclass ~= 19) = 1;
		new_peakclass(peakclass == 5) = 2;
		new_peakclass(peakclass == 14) = 3;
		new_peakclass(peakclass == 19) = 4;

						%%%TODO: You should annotate what you are doing here
						%%% 	I will have a hard time understanding if you don't annotate here
						%%% 	and this part of code somehow changed the color of my 2nd derivative plot
						%%%	also in general such modifications (any variable you add, for instance)
						%%% 	should go under File history. This way I can find out where you change
						%%%	things just by CTRL+F the variables
if find(remove_cells == -1) > 0			% above/below threshold RMSE plot
	if length(peakclass_labels) >= 10 && ...	% length over 10 should only be in cases of RMSE plots
	   length(unique(peakclass)) ~= length(peakclass_labels)	% if peakclasses are removed/absent from peakclass
		new_peakclass = reorder(new_peakclass);		% see reorder.m in /home/Matlab/Brians_Functions/ e.g. [1 2 4 4 5] -> [1 2 3 3 4]
	end
%}
