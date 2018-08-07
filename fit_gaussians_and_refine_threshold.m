function [modelBest, numComponents, muBest, stdevBest, proportionBest, minAIC, newThr1, newThr2, newPeakclass, minThreshold, gaussianOrder] = fit_gaussians_and_refine_threshold (data, fileName, label, varargin)
%% Fits data to Gaussian mixture models and finds the optimal number of components
% Usage: [modelBest, numComponents, muBest, stdevBest, proportionBest, minAIC, newThr1, newThr2, newPeakclass, minThreshold] = fit_gaussians_and_refine_threshold (data, fileName, label, varargin)
% Arguments: TODO
%       data must be a column vector of values; 
%       maxNumComponents is defaulted to be 5
%       the rest are only necessary if the histogram is to be plotted and saved
%       Examine cell plot must have figname rmse_*R/F*_row_Fit_traces and 
%           threshold plot must have figname rmse_*R/F*_row_Fit_threshold
% 
% Requires:    
%        /home/Matlab/Adams_Functions/histg.m
% Used by:    
%        /media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%        /media/adamX/m3ha/data_dclamp/find_initial_slopes.m
%
% 2016-08-03 - Created
% 2016-09-15 - Added fitmode
% 2016-09-30 - Added peakclass & peakclassLabels
% 2016-09-30 - Took out spontaneous spikes 
% 2016-11-01 - Replaced ProbDistUnivParam with makedist, also from the Statistics and Machine Learning Toolbox
% 2017-02-08 - BT - Changed legend location for RMSE graph fitting
% 2017-02-22 - BT - Identified RMSEs above threshold, returns modified peakclass for RMSE
% 2017-04-29 - Don't apply reorder.m unless plotting RMSE
% 2017-05-01 - BT - remove_cells for cells to be examined between RMSE plots. 
% 2017-05-05 - BT - Examine cell plot must have figname rmse_*R/F*_row_Fit_traces 
%               and threshold plot must have figname rmse_*R/F*_row_Fit_threshold
% 2018-01-24 - Added isdeployed
% 2018-07-18 - BT - Implemented input parser

%% Set parameters
prec = 10^-4;               % Precision
remove_cells = [5, 14, 19]; % Cell #s to examine in RMSE plotting

%% Default values for optional arguments
defaultMaxNumComponents = 5;
defaultOldThr = 0;
defaultThrMin = 0;
defaultOutFolder = pwd;
defaultFitMode = 0;
defaultPeakClass = [];
defaultPeakClass_Labels = {}; % length must be the same as unique elements in peakclass, set 'data' as prefix in dependent argument
defaultMin_Threshold = 0;
defaultThresMode = 'threeStdFirstComponent';     % 'minFirstTwoComponents' - Passive fitting for Narrowest_peak_2ndder_nospont & RMSE
                                                % 'threeStdFirstComponent' - Initial slopes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if ~isdeployed
    if exist('/home/Matlab/', 'dir') == 7
        functionsdirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsdirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsdirectory does not exist!');
    end
    addpath(fullfile(functionsdirectory, '/Adams_Functions/'));    % for histg.m
end

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'data', ...
    @(x) validateattributes(x, {'numeric'}, {'nonempty', 'column'}));
addRequired(iP, 'fileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'label', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MaxNumComponents', defaultMaxNumComponents, ...
    @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'positive', 'integer'}));
addParameter(iP, 'OldThr', defaultOldThr, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ThrMin', defaultThrMin, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'OutFolder', defaultOutFolder, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FitMode', defaultFitMode, ...
    @(x) validateattributes(x, {'numeric'}, ...
        {'>=', 0, '<=', 2}));
addParameter(iP, 'PeakClass', defaultPeakClass, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'PeakClassLabels', defaultPeakClass_Labels, ...
    @(x) validateattributes(x, {'cell'}, {'vector'}));
addParameter(iP, 'MinThreshold', defaultMin_Threshold, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ThresMode', defaultThresMode, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, data, fileName, label, varargin{:});
maxNumComponents = iP.Results.MaxNumComponents;
oldThr = iP.Results.OldThr;
thrMin = iP.Results.ThrMin;
outfolder = iP.Results.OutFolder;
fitmode = iP.Results.FitMode;
peakclass = iP.Results.PeakClass;
peakclassLabels = iP.Results.PeakClassLabels;
minThreshold = iP.Results.MinThreshold;
thresMode = iP.Results.ThresMode; 

% Set dependent argument defaults
if isempty(peakclass) && isempty(peakclassLabels)
    peakclass = ones(size(data));
    peakclassLabels = {'data'};
elseif ~isempty(peakclass) && isempty(peakclassLabels)
    peakclass_labels_prefix = peakclassLabels
    peakclassLabels = cell(size(peakclass));
    for i = 1:length(peakclass)
        peakclassLabels{i} = ['data' num2str(i)];
    end
end

%% Set up log file
logSuffix = ['_fit_gaussians_', num2str(maxNumComponents), 'maxcomplog.txt'];
logfile = fullfile(outfolder, strrep(fileName, '.png', logSuffix));
fid = fopen(logfile, 'w');

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
binw = 10^(ordomag - 1);            % Bin width
left = floor(dmin/binw)*binw;       % left bound
right = ceil(dmax/binw)*binw;       % right bound
edges = (left:binw:right)';

%% Find area of histogram
harea = binw * npts;

%% Find the best Gaussian mixture fit
if nargin < 2
    maxNumComponents = 5;
end
AIC = zeros(maxNumComponents,1);
GMModels = cell(maxNumComponents,1);   % Preallocation of Gaussian mixture models
Mus = cell(maxNumComponents,1);        % Preallocation of means of Gaussian mixture components
Stdevs = cell(maxNumComponents,1);     % Preallocation of standard deviations of Gaussian mixture components
Props = cell(maxNumComponents,1);      % Preallocation of proportions of Gaussian mixture components
options = statset('MaxIter', 1000);
rng(1);                    % For reproducibility
for k = 1:maxNumComponents
    try 
        GMModels{k} = fitgmdist(data, k, 'Options', options);    
        Mu = zeros(k, 1);
        Stdev = zeros(k, 1);
        Prop = zeros(k, 1);
        for i = 1:k
            Mu(i) = GMModels{k}.mu(i);
            Stdev(i) = sqrt(GMModels{k}.Sigma(1, 1, i));
            Prop(i) = GMModels{k}.ComponentProportion(1, i);
        end
        [Mus{k}, I_s] = sort(Mu);
        Stdevs{k} = Stdev(I_s);
        Props{k} = Prop(I_s);
        fprintf(fid, '\n GM Mean(s), Standard deviation(s) & Proportions(s) for %i Component(s)\n', k);
        for i = 1:k
            fprintf(fid, '%g\t%g\t%g\n', Mus{k}(i), Stdevs{k}(i), Props{k}(i));
        end
        AIC(k)= GMModels{k}.AIC;
    catch ME
        switch ME.identifier
            case 'stats:gmdistribution:IllCondCovIter'
                AIC(k) = NaN;
                GMModels{k} = NaN;
                Mus{k} = NaN;
                Stdevs{k} = NaN;
                Props{k} = NaN;
                warning(['Ill-conditioned covariance created, skipping ', ...
                        num2str(k), ' components.']);
            otherwise
                AIC(k) = NaN;
                GMModels{k} = NaN;
                Mus{k} = NaN;
                Stdevs{k} = NaN;
                Props{k} = NaN;
                continue;
        end
    end
end
[minAIC, numComponents] = min(AIC);
modelBest = GMModels{numComponents};
muBest = Mus{numComponents};
stdevBest = Stdevs{numComponents};
proportionBest = Props{numComponents};
fprintf(fid, '\n GM Model selected has %d components \n', numComponents);
fprintf(fid, '\n Mean(s)\tStandard deviation(s)\tProportions(s)\n');
for i = 1:numComponents
    fprintf(fid, '%g\t%g\t%g\n', muBest(i), stdevBest(i), proportionBest(i));
end

%% Extract info about the Gaussian mixture fit
x = (left:prec:right)';
Pdf = pdf(modelBest, x);
Normal = cell(numComponents, 1);
Normalpdf = cell(numComponents, 1);
for i = 1:numComponents
    % Create a normal distribution for this component
    Normal{i} = makedist('Normal', 'mu', muBest(i), 'sigma', stdevBest(i));

    % Get the probability density function
    Normalpdf{i} = pdf(Normal{i}, x);
end

%% Find new thresholds
if strcmp(thresMode, 'minFirstTwoComponents')
    if numComponents == 3                    % 2 possible thresholds

        % Find the minimum between lower two Gaussians
        left1 = max(floor(muBest(1)/prec)*prec, thrMin);  % left bound; threshold can't be lower than thrMin
        right1 = max(ceil(muBest(2)/prec)*prec, left1);     % right bound
        xx1 = (left1:prec:right1)';
        Pdf_part1 = pdf(modelBest, xx1);
        [minprob1 minprob1_ind] = min(Pdf_part1);
        newThrLower = xx1(minprob1_ind);

        % Find the minimum between upper two Gaussians
        left2 = max(floor(muBest(2)/prec)*prec, thrMin);  % left bound; threshold can't be lower than thrMin
        right2 = max(ceil(muBest(3)/prec)*prec, left2);    % right bound
        xx2 = (left2:prec:right2)';
        Pdf_part2 = pdf(modelBest, xx2);
        [minprob2 minprob2_ind] = min(Pdf_part2);
        newThrUpper = xx2(minprob2_ind);

        % Use the closest value to the old threshold as the new threshold
        [~, choice] = min(abs([newThrLower, newThrUpper] - oldThr));
        if choice == 1
            newThr1 = newThrLower;
            newThr2 = newThrUpper;
        elseif choice == 2
            newThr1 = newThrUpper;
            newThr2 = newThrLower;
        else
            error('threshold is wrong!')
        end

        fprintf(fid, 'New LTS threshold is %g V^2/s^2\n', newThr1);
        fprintf(fid, 'Alternate LTS threshold is %g V^2/s^2\n', newThr2);
        fclose(fid);
     
    % elseif numComponents == 4                % 3 possible thresholds
    else
        comp1 = find(muBest < oldThr, 1, 'last');         % Gaussian peak left of oldThr
        comp2 = find(muBest > oldThr, 1);                 % Gaussian peak right of oldThr
        value1 = muBest(comp1);
        value2 = muBest(comp2);

        left1 = max(floor(value1/prec)*prec, thrMin);      % left bound; threshold can't be lower than thrMin
        right1 = ceil(value2/prec)*prec;                    % right bound
        xx = (left1:prec:right1)';
        Pdf_part = pdf(modelBest, xx);
        [minprob minprob_ind] = min(Pdf_part);
        newThr1 = xx(minprob_ind);
        fprintf(fid, 'New LTS threshold is %g V^2/s^2\n', newThr1);
        fclose(fid);
    end
elseif strcmp(thresMode, 'threeStdFirstComponent')
    maxGaussians = zeros(size(proportionBest)); % holds maximum value per component
    for idx = 1:length(Normalpdf)   % compute maximum value of each component
        maxGaussians(idx) = max(Normalpdf{idx} * proportionBest(idx) * harea);
    end
    [~, modeIdx] = max(maxGaussians);    % index of most representative Gaussian
    newThr1 = muBest(modeIdx) + 3 * stdevBest(modeIdx);   % right bound
    newThr2 = muBest(modeIdx) - 3 * stdevBest(modeIdx);   % left bound
    gaussianOrder = zeros(3,1); % order gaussians 1) mode, 2) leftmost if mode is center, 3) rightmost if mode is center
    gaussianOrder(1) = modeIdx; % or 1) mode, 2) center if mode is leftmost, 3) rightmost if mode is leftmost
    if modeIdx == 2
        gaussianOrder(2) = 1;
        gaussianOrder(3) = 3;
    elseif modeIdx == 1
        gaussianOrder(2) = 2;
        gaussianOrder(3) = 3;
    end
end

%% Plot and save histogram with fitted Gaussian mixture pdf
if fitmode == 0
    fitsuffix = ' (all)';
elseif fitmode == 1
    fitsuffix = ' (100%, 200%, 400% g incr)';
elseif fitmode == 2
    fitsuffix = ' (for fitting)';
end

% Copy peakclass to newPeakclass as template
newPeakclass = peakclass;

% Copy peakclassLabels to new_peakclass_labels as template if not changed later
new_peakclass_labels = peakclassLabels;

if nargin > 3       %%% TODO: Why 3? BT - in the original source
    h = figure(400);
    set(h, 'Visible', 'on');
    set(h, 'Name', ['Distribution of ', label]);
    clf(h)
    opt.edges = edges;                      % needed for histg
    opt.group_names = peakclassLabels';    % needed for histg
    if length(peakclassLabels) >= 10 && ...
        strcmp(label, 'RMSE (mV) in the falling phase');
                                            % RMSE falling phase plot analysis plots
        if ~isempty(strfind(fileName, 'rmse_F_row_Fit_threshold'))        
                                            % marker for plotting above/below threshold 
            % Labels for above/below threshold plot
            new_peakclass_labels = {'Above Threshold', 'Below Threshold'};
            opt.group_names = new_peakclass_labels';        % for histg.m
        elseif ~isempty(strfind(fileName, 'rmse_F_row_Fit_traces')) && ...
                size(remove_cells,2) > 0    % plotting specific cells
%            new_peakclass_labels = {'All other cells'};        % Make labels for cells to remove
%            for q = 1:size(remove_cells,2)
%                new_peakclass_labels{q+1} = ['Cell #' num2str(remove_cells(q))];
%            end
%            opt.group_names = new_peakclass_labels';       % for histg.m
            % Make labels for cells to remove
            new_peakclass_labels = {'All other cells'};
            not_remove_cells_indices = {};                  % indices of cells not in remove_cells
            for q = 1:size(remove_cells,2)
                new_peakclass_labels{q+1} = ['Cell #' num2str(remove_cells(q))];    % make cell labels
                not_remove_cells_indices{q} = find(peakclass ~= remove_cells(q));    
                                % track indices of cells not being examined
            end
            opt.group_names = new_peakclass_labels';        % for histg.m
            newPeakclass(intersectm(not_remove_cells_indices)) = 1;    
                                % change all indices of cells not being examined to 1
            for q = 1:size(remove_cells,2)  % cells being examined given new classes for histg.m
                newPeakclass(peakclass == remove_cells(q)) = q + 1;
            end
        end
    end
    z = histg(data, newPeakclass, opt);    % plots stacked histogram
    for i = 1:numComponents                 % plots gaussian components
        plot(x, Normalpdf{i} * proportionBest(i) * harea, 'w', 'Displayname', ['component #', num2str(i)]);
    end
    full_dist = plot(x, Pdf * harea, 'r', 'LineWidth', 1, ...
                    'Displayname', 'full distribution');    % plots full distribution
    if numComponents == 3
        line([newThr1 newThr1], [0 npts], 'Color', 'r', 'LineStyle', '--', ...
                'Displayname', ['New threshold = ', num2str(newThr1)]); hold on;
        line([newThr2 newThr2], [0 npts], 'Color', 'b', 'LineStyle', '--', ...
            'Displayname', ['Alternate threshold = ', num2str(newThr2)]);
    end
    % newThr1 = -1; newThr2 = -1;              % temp        %%% TODO: annotate here
    % xlim([left right]); 
    ylim([0 max(Pdf * harea)]); 
    xlabel(label)
    ylabel('# of sweeps')
    title(['Distribution of ', label, fitsuffix]);
    if length(peakclassLabels) < 10        % Legend for 2nd derivative, using peakclass
        % legend('Location', 'northwest'); %
        legend('Location', 'northeast'); %%% for fitting initial slopes
    elseif length(peakclassLabels) >= 10 && ...
        strcmp(label, 'RMSE (mV) in the rising phase')
        % Legend for rising RMSE, location northwest overlaps plot; peakclass = cellidrow
        [min_val, min_ind] = min(full_dist.YData);  % Minimum y value of data to determine threshold
        x_min = full_dist.XData(min_ind);           % Corresponding x value of minimum y value
        if strfind(fileName, 'rmse_R_row_Fit_threshold')    % above/below threshold RMSE plot
            newPeakclass(data >= x_min) = 2;       % Change peakclasses to visibly identify traces above threshold
            newPeakclass(data < x_min) = 1;        % All data above and below min x value altered
            new_peakclass_labels = {'Above Threshold', 'Below Threshold'};
            opt.group_names = new_peakclass_labels';
        elseif ~isempty(strfind(fileName, 'rmse_R_row_Fit_traces')) && ...
                size(remove_cells,2) > 0        % cells are being examined between rising/falling plots
            new_peakclass_labels = {'All other cells'};        % Make labels for cells to remove
            not_remove_cells_indices = {};                % indices of cells not in remove_cells
            for q = 1:size(remove_cells,2)
                new_peakclass_labels{q+1} = ['Cell #' num2str(remove_cells(q))];    % make cell labels
                not_remove_cells_indices{q} = find(peakclass ~= remove_cells(q));    
                                % track indices of cells not being examined
            end
            opt.group_names = new_peakclass_labels';    % for histg.m
            newPeakclass(intersectm(not_remove_cells_indices)) = 1;    
                                % change all indices of cells not being examined to 1
            for q = 1:size(remove_cells,2)        % cells being examined given new classes for histg.m
                newPeakclass(peakclass == remove_cells(q)) = q+1;
            end
        end
        clf(h);    % Clear figure to replot with new peakclasses
        histg(data, newPeakclass, opt);    % Replot with new peakclass identifier
        newPeakclass = newPeakclass';        % Inverted to pass as argument from rising to falling phase
        for i = 1:numComponents            % Plot each Gaussian component
            plot(x, Normalpdf{i} * proportionBest(i) * harea, 'w', 'Displayname', ['component #', num2str(i)]);
        end
        full_dist = plot(x, Pdf * harea, 'r', 'LineWidth', 1, 'Displayname', 'full distribution');    % Plot the final fit
        xlim([left right]); 
        ylim([0 max(Pdf * harea)]); 
        xlabel(label)
        ylabel('# of sweeps')
        title(['Distribution of ', label, fitsuffix]);
        legend(new_peakclass_labels, 'Location', 'eastoutside');    
                                        % Revised legend, original if no analysis criteria met
%        legend('Location', 'eastoutside');    % Revised legend, original if no analysis criteria met
        line([x_min x_min], [0 npts], 'Color', 'r', 'LineStyle', '--', ...
            'Displayname', ['New threshold = ', num2str(newThr1)]); hold on;    % Plot threshold line at minimum in rising
        text(1.1, 0.95, ['Rising Threshold: ' num2str(x_min)], 'Units', 'normalized');        
                                        % Text at what RMSE the threshold is at
            text(1.1, 0.85, ['Components: ' num2str(maxNumComponents)], 'Units', 'normalized');    
                                        % Text how many components are present
        minThreshold = x_min;    % Pass as argument to falling phase for labelling
      elseif length(peakclassLabels) >= 10 && ...
        strcmp(label, 'RMSE (mV) in the falling phase')        % Legend for falling RMSE without replotting
        legend(new_peakclass_labels, 'Location', 'eastoutside');    % Revised legend, original if no analysis criteria met
%        legend('Location', 'eastoutside');    % Revised legend, original if no analysis criteria met
        text(1.1, 0.95, ['Rising Threshold: ' num2str(minThreshold)], 'Units', 'normalized');        % Text at what RMSE the threshold is at
        text(1.1, 0.85, ['Components: ' num2str(maxNumComponents)], 'Units', 'normalized');    % Text how many components are present    
    end
end
if nargin > 4        %%% TODO: Why 4? BT - also in original source
    figname = fullfile(outfolder, fileName);
    saveas(h, figname);
    xlim([-0.009 0]);
    ylim([0 500]);
    title(['Distribution of ', label, ', zoomed in', fitsuffix]);
    figname = fullfile(outfolder, strrep(fileName, '.png', '_zoomed.png'));
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
ind_cl = cell(1, 7);                    % Indices for each peak class
data_cl = cell(1, 7);                    % data for each peak class
color_cl = {'y', 'm', 'c', 'r', 'g', 'b', 'k'};
for pkcl = 1:7
    ind_cl{pkcl} = find(peakclass == pkcl);
    data_cl{pkcl} = data(ind_cl{pkcl});
end

    for pkcl = 1:7
        histogram(data_cl{pkcl}, edges, 'FaceColor', color_cl{pkcl}, 'EdgeColor', color_cl{pkcl});
    end

    if numComponents == 3
        legend(['New threshold = ', num2str(newThr1)], ...
            ['Alternate threshold = ', num2str(newThr2)], ...
            'Peak class #1', 'Peak class #2', 'Peak class #3', ...
            'Peak class #4', 'Peak class #5', 'Peak class #6', ...
            'Peak class #7', 'Location', 'northwest');
    else
        legend(['New threshold = ', num2str(newThr1)], ...
            'Peak class #1', 'Peak class #2', 'Peak class #3', ...
            'Peak class #4', 'Peak class #5', 'Peak class #6', ...
            'Peak class #7', 'Location', 'northwest');
    end

Normal{i} = ProbDistUnivParam('normal', [muBest(i) stdevBest(i)]);
newPeakclass(peakclass ~= 5 & peakclass ~= 14 & peakclass ~= 19) = 1;
        newPeakclass(peakclass == 5) = 2;
        newPeakclass(peakclass == 14) = 3;
        newPeakclass(peakclass == 19) = 4;

                        %%%TODO: You should annotate what you are doing here
                        %%%     I will have a hard time understanding if you don't annotate here
                        %%%     and this part of code somehow changed the color of my 2nd derivative plot
                        %%%    also in general such modifications (any variable you add, for instance)
                        %%%     should go under File history. This way I can find out where you change
                        %%%    things just by CTRL+F the variables
if find(remove_cells == -1) > 0            % above/below threshold RMSE plot
    if length(peakclassLabels) >= 10 && ...    % length over 10 should only be in cases of RMSE plots
       length(unique(peakclass)) ~= length(peakclassLabels)    % if peakclasses are removed/absent from peakclass
        newPeakclass = reorder(newPeakclass);        % see reorder.m in /home/Matlab/Brians_Functions/ e.g. [1 2 4 4 5] -> [1 2 3 3 4]
    end
addpath(fullfile(functionsdirectory, '/Brians_Functions/'));    % for reorder.m
%        /home/Matlab/Brians_Functions/reorder.m
function [modelBest, numComponents, muBest, stdevBest, proportionBest, minAIC, newThr1, newThr2, newPeakclass, minThreshold] = fit_gaussians_and_refine_threshold (data, maxNumComponents, oldThr, thrMin, label, outfolder, fileName, fitmode, peakclass, peakclassLabels, minThreshold)

data = iP.Results.Data;
fileName = iP.Results.FileName;
label = iP.Results.Label;

%}
