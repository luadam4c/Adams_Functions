function fits = compute_iei_thresholds (data, outFolderParent)
%% Compute inter-event-interval thresholds that separates events from spikes
% 
% Requires:
%       cd/nanstderr.m
%       cd/color_index.m
%       cd/fit_logIEI.m
%       cd/bar_w_CI.m
%       cd/plot_grouped_histogram.m
%       cd/plot_grouped_scatter.m
%
% File History:
% 2018-07-25 Adapted code from zgPeriodStats.m, 
%               which was derived from paula_iei4.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
xUnitHist = 'log s';
plotFitsFlag = true;        % whether to plot individual fits
plotOverlapFlag = false;    % whether to plot overlapped fits across slices/cells
plotHeatMapFlag = false;    % whether to plot fits across slices/cells as a heatmap
plotBarFlag = false;        % whether to plot barplots of fitted parameters
plotStatHistFlag = false;   % whether to plot histograms of fitted parameters
plotScatterFlag = false;    % whether to plot scatter plots
plotPerCellFlag = false;    % whether to plot time-consuming plots per cell
truncateDataFlag = false;   % whether to truncate data
nPoints = 500;              % number of points for plotting pdfs
maxPdfValue = 4.5;          % maximum probability density value
cutPdfValue = 0.1;          % cutoff probability density value

%% Set the order of experiments
allExpts = {'xHC067_3nMAng_Control', ...
            'x10HC067_3nMAng', ...
            'x50pMAng', 'x300pMAng', ...
            'x3nMAng', 'x1uMang'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create outFolderParent if it doesn't exist
if exist(outFolderParent, 'dir') ~= 2
    mkdir(outFolderParent);
end

%% Import and extract logIEIs
perExptIEIs = data.perExpt.periods;

%%%% Per Experiment

%% Fit and plot Paula's IEI logIEIs
% Determine the experiments to fit
allExpts = fieldnames(perExptIEIs);     % Replace by defined order above
nExpts = numel(allExpts);

% Fit and plot pooled logIEIs for each experiment
muModel0 = zeros(2, nExpts);
muModel2 = zeros(2, nExpts);
muModel2Low = zeros(2, nExpts);
muModel2High = zeros(2, nExpts);
muModel3 = zeros(2, nExpts);
muModel3Low = zeros(2, nExpts);
muModel3High = zeros(2, nExpts);
minLogIEIs = zeros(nExpts, 1);
maxLogIEIs = zeros(nExpts, 1);
for iExpt = 1:nExpts            % for each experiment
    % Extract identifier for experiment
    idExpt = allExpts{iExpt};
    fprintf('Experiment %s ... \n', idExpt);

    % Extract all IEIs pooled from different slices of this experiment
    IEIs = perExptIEIs.(idExpt).allCatted;
    IEIs = IEIs(IEIs > 0);
    logIEIs = log(IEIs);


    % Determine maximum and minimum IEIs
    if ~isempty(logIEIs)
        minLogIEIs(iExpt) = min(logIEIs);
        maxLogIEIs(iExpt) = max(logIEIs);
    else
        minLogIEIs(iExpt) = NaN;
        maxLogIEIs(iExpt) = NaN; 
    end
 
    fprintf('logIEIs extracted from %s ... \n', idExpt);

    % Fit logIEIs to distributions and plot on top of histograms
    if ~isempty(logIEIs)
        [~, ~, ~, ~, params] = ...
                fit_logIEI(logIEIs, idExpt, 'XUnit', xUnitHist, ...
                        'TruncateFlag', truncateDataFlag, ...
                        'PlotFlag', plotFitsFlag, ...
                        'OutFolder', outFolderParent);
        fprintf('Done with experiment %s! \n\n', idExpt);

        % Store parameter fits in data structure
        fits.(idExpt) = params;

        % Extract information from fits
        if isfield(params, 'mu2Model0')
            muModel0(1, iExpt) = exp(params.mu1Model0);
            muModel0(2, iExpt) = exp(params.mu2Model0);
        else
            muModel0(1, iExpt) = NaN;
            muModel0(2, iExpt) = NaN;
        end
        if isfield(params, 'mu1Model2')
            muModel2(1, iExpt) = exp(params.mu1Model2);
            muModel2Low(1, iExpt) = exp(params.mu1CIModel2(1));
            muModel2High(1, iExpt) = exp(params.mu1CIModel2(2));
            muModel2(2, iExpt) = exp(params.mu2Model2);
            muModel2Low(2, iExpt) = exp(params.mu2CIModel2(1));
            muModel2High(2, iExpt) = exp(params.mu2CIModel2(2));
        else
            muModel2(1, iExpt) = NaN;
            muModel2Low(1, iExpt) = NaN;
            muModel2High(1, iExpt) = NaN;
            muModel2(2, iExpt) = NaN;
            muModel2Low(2, iExpt) = NaN;
            muModel2High(2, iExpt) = NaN;
        end
        if isfield(params, 'mu1Model3')
            muModel3(1, iExpt) = exp(params.mu1Model3);
            muModel3Low(1, iExpt) = exp(params.mu1CIModel3(1));
            muModel3High(1, iExpt) = exp(params.mu1CIModel3(2));
            muModel3(2, iExpt) = params.mu2Model3;
            muModel3Low(2, iExpt) = params.mu2CIModel3(1);
            muModel3High(2, iExpt) = params.mu2CIModel3(2);
        else
            muModel3(1, iExpt) = NaN;
            muModel3Low(1, iExpt) = NaN;
            muModel3High(1, iExpt) = NaN;
            muModel3(2, iExpt) = NaN;
            muModel3Low(2, iExpt) = NaN;
            muModel3High(2, iExpt) = NaN;
        end
    else
        fprintf('No logIEIs for experiment %s!\n\n', idExpt);
    end
end

% Compute x limits for plots; each row is an experiment
xlimits = [minLogIEIs, maxLogIEIs];

% Compute x vectors for plots; each row is an experiment
xVecs = cell(nExpts, 1);
parfor iExpt = 1:nExpts
    xVecs{iExpt} = linspace(xlimits(iExpt, 1), xlimits(iExpt, 2), nPoints);
end

% Plot bar graphs comparing means
if plotBarFlag
    h = figure(1004);
    clf(h);
    bar(muModel0);
    set(gca, 'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using a kernel distribution for pooled logIEIs');
    saveas(h, fullfile(outFolderParent, 'Kernel_means_Pooled'), 'png');

    h = figure(1005);
    clf(h);
    h = bar_w_CI(h, muModel2, muModel2Low, muModel2High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using two-Gaussian fits for pooled logIEIs');
    saveas(h, fullfile(outFolderParent, 'GaussOnly_means_Pooled'), 'png');

    h = figure(1006);
    clf(h);
    h = bar_w_CI(h, muModel3, muModel3Low, muModel3High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using Gaussian-Exp-Exponential fits for pooled logIEIs');
    saveas(h, fullfile(outFolderParent, 'GaussExpExp_means_Pooled'), 'png');

    h = figure(1007);
    clf(h);
    h = bar_w_CI(h, [muModel2(1, :); muModel3(1, :)], ...
                    [muModel2Low(1, :); muModel3Low(1, :)], ...
                    [muModel2High(1, :); muModel3High(1, :)], ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Mean inter-spike interval (s)');
    title('Mean intraburst ISIs for pooled logIEIs');
    saveas(h, fullfile(outFolderParent, 'Mu1_Pooled'), 'png');
end

%%%% Per Slice

% Fit and plot logIEIs from each slice for each experiment
muPerSliceModel0 = zeros(2, nExpts);    % means averaged from each slice
muPerSliceModel0Low = zeros(2, nExpts);
muPerSliceModel0High = zeros(2, nExpts);
muPerSliceModel2 = zeros(2, nExpts);    % means averaged from each slice
muPerSliceModel2Low = zeros(2, nExpts);
muPerSliceModel2High = zeros(2, nExpts);
muPerSliceModel3 = zeros(2, nExpts);    % means averaged from each slice
muPerSliceModel3Low = zeros(2, nExpts);
muPerSliceModel3High = zeros(2, nExpts);
spacingPerSlice = zeros(2, nExpts);     % spacing parameters averaged from each slice
spacingPerSliceLow = zeros(2, nExpts);
spacingPerSliceHigh = zeros(2, nExpts);
threshold1PerSlice = zeros(2, nExpts);  % threshold #1s averaged from each slice
threshold1PerSliceLow = zeros(2, nExpts);
threshold1PerSliceHigh = zeros(2, nExpts);
threshold2PerSlice = zeros(2, nExpts);  % threshold #2s averaged from each slice
threshold2PerSliceLow = zeros(2, nExpts);
threshold2PerSliceHigh = zeros(2, nExpts);
void1PerSlice = zeros(2, nExpts);       % void parameter #1s averaged from each slice
void1PerSliceLow = zeros(2, nExpts);
void1PerSliceHigh = zeros(2, nExpts);
void2PerSlice = zeros(2, nExpts);       % void parameter #2s averaged from each slice
void2PerSliceLow = zeros(2, nExpts);
void2PerSliceHigh = zeros(2, nExpts);
for iExpt = 1:nExpts            % for each experiment
    % Extract identifier for experiment
    idExpt = allExpts{iExpt};
    if plotOverlapFlag || plotHeatMapFlag
        xVec = xVecs{iExpt};
    end
    fprintf('Experiment %s ... \n', idExpt);

    % Extract IEIs for this experiment and remove pooled logIEIs
    thisExpt = rmfield(perExptIEIs.(idExpt), 'allCatted');  % a structure

    % Get all slices
    idAllSlices = fieldnames(thisExpt);
    nSlices = numel(idAllSlices);
    ctSliceModel0 = 0;                  % counts useful slices for Model 0
    ctSliceModel2 = 0;                  % counts useful slices for Model 2
    ctSliceModel3 = 0;                  % counts useful slices for Model 3
    mu1Model0AllSlices = [];            % stores means of intra-burst intervals
    mu2Model0AllSlices = [];            % stores means of inter-burst intervals
    spacingModel0AllSlices = [];        % stores spacing parameters
    threshold1Model0AllSlices = [];     % stores burst detection threshold #1s
    threshold2Model0AllSlices = [];     % stores burst detection threshold #2s
    void1Model0AllSlices = [];          % stores void parameter #1s
    void2Model0AllSlices = [];          % stores void parameter #2s
    mu1Model2AllSlices = [];            % stores means of intra-burst intervals
    mu2Model2AllSlices = [];            % stores means of inter-burst intervals
    spacingModel2AllSlices = [];        % stores spacing parameters
    threshold1Model2AllSlices = [];     % stores burst detection threshold #1s
    threshold2Model2AllSlices = [];     % stores burst detection threshold #2s
    void1Model2AllSlices = [];          % stores void parameter #1s
    void2Model2AllSlices = [];          % stores void parameter #2s
    mu1Model3AllSlices = [];            % stores means of intra-burst intervals
    mu2Model3AllSlices = [];            % stores means of inter-burst intervals
    spacingModel3AllSlices = [];        % stores spacing parameters
    threshold1Model3AllSlices = [];     % stores burst detection threshold #1s
    threshold2Model3AllSlices = [];     % stores burst detection threshold #2s
    void1Model3AllSlices = [];          % stores void parameter #1s
    void2Model3AllSlices = [];          % stores void parameter #2s
    if plotOverlapFlag || plotHeatMapFlag
        pdfModel0AllSlices = [];        % stores pdfs for Model 0
        pdfModel2AllSlices = [];        % stores pdfs for Model 2
        comp2Model2AllSlices = [];      % stores component #2 pdfs for Model 2
        pdfModel3AllSlices = [];        % stores pdfs for Model 3
        comp2Model3AllSlices = [];      % stores component #2 pdfs for Model 3
    end

    % Construct outFolder for this experiment
    outFolder = fullfile(outFolderParent, idExpt);
    if exist(outFolder, 'dir') ~= 7 && ...
        (plotBarFlag  || plotFitsFlag || plotOverlapFlag || plotHeatMapFlag)
        mkdir(outFolder);
    end

    for iSlice = 1:nSlices      % for each slice
        % Extract identifier for slice
        idSlice = idAllSlices{iSlice};
        idSliceFull = [idExpt, '_', idSlice];
        fprintf('Slice %s ... \n', idSlice);

        % Extract IEIs pooled from all cells of this slice
        IEIs = thisExpt.(idSlice).allCells;
        IEIs = IEIs(IEIs > 0);
        logIEIs = log(IEIs);
        fprintf('logIEIs extracted from %s ... \n', idSliceFull);

        % Fit logIEIs to distributions and plot on top of histograms
        if ~isempty(logIEIs)
            [~, ~, ~, ~, params] = ...
                    fit_logIEI(logIEIs, idSliceFull, 'XUnit', xUnitHist, ...
                            'TruncateFlag', truncateDataFlag, ...
                            'PlotFlag', plotFitsFlag, ...
                            'OutFolder', outFolder);

            % Store parameter fits in data structure
            fits.(idExpt).(idSlice) = params;

            if isfield(params, 'mu2Model0')
                ctSliceModel0 = ctSliceModel0 + 1;
                mu1Model0AllSlices(ctSliceModel0, 1) = params.mu1Model0;
                mu2Model0AllSlices(ctSliceModel0, 1) = params.mu2Model0;
                spacingModel0AllSlices(ctSliceModel0, 1) = params.spacingModel0;
                threshold1Model0AllSlices(ctSliceModel0, 1) = params.thresholdModel0;
                threshold2Model0AllSlices(ctSliceModel0, 1) = params.thresholdModel0;
                void1Model0AllSlices(ctSliceModel0, 1) = params.voidModel0;
                void2Model0AllSlices(ctSliceModel0, 1) = params.voidModel0;
                if plotOverlapFlag || plotHeatMapFlag
                    pdfModel0AllSlices = [pdfModel0AllSlices; ...
                                            params.pdfModel0(xVec)];
                end
                fprintf('Done fitting %s with a kernel distribution! \n\n', idSliceFull);
            else
                fprintf('%s was not fitted with a kernel distribution! \n\n', idSliceFull);
            end
            if isfield(params, 'mu1Model2')
                ctSliceModel2 = ctSliceModel2 + 1;
                mu1Model2AllSlices(ctSliceModel2, 1) = params.mu1Model2;
                mu2Model2AllSlices(ctSliceModel2, 1) = params.mu2Model2;
                spacingModel2AllSlices(ctSliceModel2, 1) = params.spacingModel2;
                threshold1Model2AllSlices(ctSliceModel2, 1) = params.threshold1Model2;
                threshold2Model2AllSlices(ctSliceModel2, 1) = params.threshold2Model2;
                void1Model2AllSlices(ctSliceModel2, 1) = params.void1Model2;
                void2Model2AllSlices(ctSliceModel2, 1) = params.void2Model2;
                if plotOverlapFlag || plotHeatMapFlag
                    pdfModel2AllSlices = [pdfModel2AllSlices; ...
                                            params.pdfModel2(xVec)];
                    comp2Model2AllSlices = [comp2Model2AllSlices; ...
                                            params.comp2Model2(xVec)];
                end
                fprintf('Done fitting %s with two Gaussians! \n\n', idSliceFull);
            else
                fprintf('%s was not fitted with two Gaussians! \n\n', idSliceFull);
            end
            if isfield(params, 'mu1Model3')
                ctSliceModel3 = ctSliceModel3 + 1;
                mu1Model3AllSlices(ctSliceModel3, 1) = params.mu1Model3;
                mu2Model3AllSlices(ctSliceModel3, 1) = params.mu2Model3;
                spacingModel3AllSlices(ctSliceModel3, 1) = params.spacingModel3;
                threshold1Model3AllSlices(ctSliceModel3, 1) = params.threshold1Model3;
                threshold2Model3AllSlices(ctSliceModel3, 1) = params.threshold2Model3;
                void1Model3AllSlices(ctSliceModel3, 1) = params.void1Model3;
                void2Model3AllSlices(ctSliceModel3, 1) = params.void2Model3;
                if plotOverlapFlag || plotHeatMapFlag
                    pdfModel3AllSlices = [pdfModel3AllSlices; ...
                                            params.pdfModel3(xVec)];
                    comp2Model3AllSlices = [comp2Model3AllSlices; ...
                                            params.comp2Model3(xVec)];
                end
                fprintf('Done fitting %s with Gaussian-Exp-Exponential! \n\n', ...
                        idSliceFull);
            else
                fprintf('%s was not fitted with Gaussian-Exp-Exponential! \n\n', ...
                        idSliceFull);
            end
        else
            fprintf('No logIEIs for %s!\n\n', idSliceFull);
        end        
    end

    % Extract information from fits
    commands = generate_stats_commands('Slice', 'Model0');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Slice', 'Model2');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Slice', 'Model3');
    cellfun(@eval, commands);

    % Plot pdfs of all slices together
    if plotOverlapFlag
        if ~isempty(pdfModel0AllSlices)
            h = figure('Visible', 'off');
            hold on;
            plot(xVec, pdfModel0AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Kernel distributions for slices of ', idExpt], ...
                    'Interpreter', 'none');
            legend(idAllSlices, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(outFolderParent, [idExpt, '_pdfs_Kernel']), 'png');
            close(h);
        end

        if ~isempty(pdfModel2AllSlices)
            h = figure('Visible', 'off');
            hold on;
            plot(xVec, pdfModel2AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Two-Gaussian fits for slices of ', idExpt], ...
                    'Interpreter', 'none');
            legend(idAllSlices, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(outFolderParent, [idExpt, '_pdfs_GaussOnly']), 'png');
            close(h);

            h = figure('Visible', 'off');
            hold on;
            plot(xVec, comp2Model2AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['2nd Gaussian fits for slices of ', idExpt], ...
                    'Interpreter', 'none');
            legend(idAllSlices, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(outFolderParent, [idExpt, '_comp2s_GaussOnly']), 'png');
            close(h);
        end

        if ~isempty(pdfModel3AllSlices)
            h = figure('Visible', 'off');
            hold on;
            plot(xVec, pdfModel3AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Gaussian-Exp-Exponential fits for slices of ', idExpt], ...
                    'Interpreter', 'none');
            legend(idAllSlices, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(outFolderParent, [idExpt, '_pdfs_GaussExpExp']), 'png');
            close(h);

            h = figure('Visible', 'off');
            hold on;
            plot(xVec, comp2Model3AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Exp-Exponential fits for slices of ', idExpt], ...
                    'Interpreter', 'none');
            legend(idAllSlices, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(outFolderParent, [idExpt, '_comp2s_GaussExpExp']), 'png');
            close(h);
        end
    end
    if plotHeatMapFlag
        commands = generate_heatmap_commands('Slice', 'Model0', outFolderParent);
        cellfun(@eval, commands);
        commands = generate_heatmap_commands('Slice', 'Model2', outFolderParent);
        cellfun(@eval, commands);
        commands = generate_heatmap_commands('Slice', 'Model3', outFolderParent);
        cellfun(@eval, commands);
    end
end

% Plot bar graph comparing means
if plotBarFlag
    h = figure(1007);
    clf(h);
    h = bar_w_CI(h, muPerSliceModel2, muPerSliceModel2Low, ...
                    muPerSliceModel2High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using two Gaussians for each slice');
    saveas(h, fullfile(outFolderParent, 'GaussOnly_means_AvgSlices'), 'png');

    h = figure(1008);
    clf(h);
    h = bar_w_CI(h, muPerSliceModel3, muPerSliceModel3Low, ...
                    muPerSliceModel3High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using Gaussian-Exp-Exponential fits for each slice');
    saveas(h, fullfile(outFolderParent, 'GaussExpExp_means_AvgSlices'), 'png');

    h = figure(1009);
    clf(h);
    h = bar_w_CI(h, [muPerSliceModel2(1, :); muPerSliceModel3(1, :)], ...
                    [muPerSliceModel2Low(1, :); muPerSliceModel3Low(1, :)], ...
                    [muPerSliceModel2High(1, :); muPerSliceModel3High(1, :)], ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Mean inter-spike interval (s)');
    title('Mean intraburst ISIs averaged across slices');
    saveas(h, fullfile(outFolderParent, 'mu1_AvgSlices'), 'png');

    h = figure(1010);
    clf(h);
    h = bar_w_CI(h, spacingPerSlice, spacingPerSliceLow, spacingPerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Time scale spacing parameter');
    title('Spacing parameters averaged across slices');
    saveas(h, fullfile(outFolderParent, 'Spacing_AvgSlices'), 'png');

    h = figure(1011);
    clf(h);
    h = bar_w_CI(h, threshold1PerSlice, threshold1PerSliceLow, threshold1PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #1 (s)');
    title('Burst detection threshold #1 averaged across slices');
    saveas(h, fullfile(outFolderParent, 'Threshold1_AvgSlices'), 'png');

    h = figure(1012);
    clf(h);
    h = bar_w_CI(h, threshold2PerSlice, threshold2PerSliceLow, threshold2PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #2 (s)');
    title('Burst detection threshold #2 averaged across slices');
    saveas(h, fullfile(outFolderParent, 'Threshold2_AvgSlices'), 'png');

    h = figure(1013);
    clf(h);
    h = bar_w_CI(h, void1PerSlice, void1PerSliceLow, void1PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #1');
    title('Void parameter #1s averaged across slices');
    saveas(h, fullfile(outFolderParent, 'Void1_AvgSlices'), 'png');

    h = figure(1014);
    clf(h);
    h = bar_w_CI(h, void2PerSlice, void2PerSliceLow, void2PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #2');
    title('Void parameter #2s averaged across slices');
    saveas(h, fullfile(outFolderParent, 'Void2_AvgSlices'), 'png');
end

%%%% Per Cell

% Fit and plot logIEIs from each cell for each experiment
muPerCellModel0 = zeros(2, nExpts);     % means averaged from each cell
muPerCellModel0Low = zeros(2, nExpts);
muPerCellModel0High = zeros(2, nExpts);
muPerCellModel2 = zeros(2, nExpts);     % means averaged from each cell
muPerCellModel2Low = zeros(2, nExpts);
muPerCellModel2High = zeros(2, nExpts);
muPerCellModel3 = zeros(2, nExpts);     % means averaged from each cell
muPerCellModel3Low = zeros(2, nExpts);
muPerCellModel3High = zeros(2, nExpts);
spacingPerCell = zeros(2, nExpts);      % spacing parameters averaged from each cell
spacingPerCellLow = zeros(2, nExpts);
spacingPerCellHigh = zeros(2, nExpts);
threshold1PerCell = zeros(2, nExpts);    % threshold #1s averaged from each cell
threshold1PerCellLow = zeros(2, nExpts);
threshold1PerCellHigh = zeros(2, nExpts);
void1PerCell = zeros(2, nExpts);         % void parameter #1s averaged from each cell
void1PerCellLow = zeros(2, nExpts);
void1PerCellHigh = zeros(2, nExpts);
for iExpt = 1:nExpts            % for each experiment
    % Extract identifier for experiment
    idExpt = allExpts{iExpt};
    if plotOverlapFlag || plotHeatMapFlag
        xVec = xVecs{iExpt};
    end
    fprintf('Experiment %s ... \n', idExpt);

    % Extract IEIs for this experiment and remove pooled logIEIs
    thisExpt = rmfield(perExptIEIs.(idExpt), 'allCatted');  % a structure

    % Get all slices
    idAllSlices = fieldnames(thisExpt);
    nSlices = numel(idAllSlices);
    ctCellModel0 = 0;               % counts useful cells
    sliceNumModel0 = [];            % stores slice number of all cells
    mu1Model0AllCells = [];         % stores means of intra-burst intervals
    mu2Model0AllCells = [];         % stores means of inter-burst intervals
    spacingModel0AllCells = [];     % stores spacing parameters
    threshold1Model0AllCells = [];  % stores burst detection threshold #1s
    threshold2Model0AllCells = [];  % stores burst detection threshold #2s
    void1Model0AllCells = [];       % stores void parameter #1s
    void2Model0AllCells = [];       % stores void parameter #2s
    ctCellModel2 = 0;               % counts useful cells
    sliceNumModel2 = [];            % stores slice number of all cells
    mu1Model2AllCells = [];         % stores means of intra-burst intervals
    mu2Model2AllCells = [];         % stores means of inter-burst intervals
    spacingModel2AllCells = [];     % stores spacing parameters
    threshold1Model2AllCells = [];  % stores burst detection threshold #1s
    threshold2Model2AllCells = [];  % stores burst detection threshold #2s
    void1Model2AllCells = [];       % stores void parameter #1s
    void2Model2AllCells = [];       % stores void parameter #2s
    ctCellModel3 = 0;               % counts useful cells
    sliceNumModel3 = [];            % stores slice number of all cells
    mu1Model3AllCells = [];         % stores means of intra-burst intervals
    mu2Model3AllCells = [];         % stores means of inter-burst intervals
    spacingModel3AllCells = [];     % stores spacing parameters
    threshold1Model3AllCells = [];  % stores burst detection threshold #1s
    threshold2Model3AllCells = [];  % stores burst detection threshold #2s
    void1Model3AllCells = [];       % stores void parameter #1s
    void2Model3AllCells = [];       % stores void parameter #2s
    if plotOverlapFlag || plotHeatMapFlag
        pdfModel0AllCells = cell(1, nSlices);     % stores pdfs for Model 0
        pdfModel2AllCells = cell(1, nSlices);     % stores pdfs for Model 2
        comp2Model2AllCells = cell(1, nSlices);   % stores component #2 pdfs for Model 2
        pdfModel3AllCells = cell(1, nSlices);     % stores pdfs for Model 3
        comp2Model3AllCells = cell(1, nSlices);   % stores component #2 pdfs for Model 3

        idAllCells = cell(nSlices); % stores IDs of all cells in the slice
    end
    for iSlice = 1:nSlices      % for each slice
        % Extract identifier for slice
        idSlice = idAllSlices{iSlice};
        fprintf('Slice %s ... \n', idSlice);

        % Construct outFolder for this slice
        outFolder = fullfile(outFolderParent, idExpt, idSlice);
        if exist(outFolder, 'dir') ~= 7 && (plotBarFlag || plotFitsFlag)
            mkdir(outFolder);
        end

        % Extract IEIs for each cell of this slice and remove pooled logIEIs
        thisSlice = thisExpt.(idSlice).perCell;             % a cell array
        nCells = numel(thisSlice);
        
        for iCell = 1:nCells    % for each cell
            % Construct identifier for this cell
            idCell = ['Cell', num2str(iCell)];
            idAllCells{iSlice}{iCell} = idCell;
            idCellFull = [idExpt, '_', idSlice, '_', idCell];
            fprintf('Cell #%d ... \n', iCell);
           
            % Extract IEIs for this cell
            IEIs = thisSlice{iCell};
            IEIs = IEIs(IEIs > 0);
            logIEIs = log(IEIs);
            fprintf('logIEIs extracted from %s ... \n', idCellFull);
  
            
% Mark's Note: I was debugging and noticed that it crashed on line 627 when there was only 1 logIEI (e.g.: x10HC067_3nMAng_data459_Cell12).
% I amended the code so that it skips logIEIs that are EITHER empty OR only
% contain one logIEI (see line 627)
            if strcmp('x10HC067_3nMAng_data459_Cell12', idCellFull) == 1
                m = 1
            end

            % Fit logIEIs to distributions and plot on top of histograms
            if ~isempty(logIEIs) && size(logIEIs,1) > 1
                [~, ~, ~, ~, params] = ...
                        fit_logIEI(logIEIs, idCellFull, 'XUnit', xUnitHist, ...
                                'TruncateFlag', truncateDataFlag, ...
                                'PlotFlag', all([plotFitsFlag, plotPerCellFlag]), ...
                                'OutFolder', outFolder);

                % Store parameter fits in data structure
                fits.(idExpt).(idSlice).(idCell) = params;

                if isfield(params, 'mu2Model0')
                    ctCellModel0 = ctCellModel0 + 1;
                    sliceNumModel0(ctCellModel0, 1) = iSlice;
                    mu1Model0AllCells(ctCellModel0, 1) = params.mu1Model0;
                    mu2Model0AllCells(ctCellModel0, 1) = params.mu2Model0;
                    spacingModel0AllCells(ctCellModel0, 1) = params.spacingModel0;
                    threshold1Model0AllCells(ctCellModel0, 1) = params.thresholdModel0;
                    threshold2Model0AllCells(ctCellModel0, 1) = params.thresholdModel0;
                    void1Model0AllCells(ctCellModel0, 1) = params.voidModel0;
                    void2Model0AllCells(ctCellModel0, 1) = params.voidModel0;
                    if plotOverlapFlag || plotHeatMapFlag
                        pdfModel0AllCells{iSlice} = [pdfModel0AllCells{iSlice}; ...
                                                        params.pdfModel0(xVec)];
                    end
                    fprintf('Done fitting %s with a kernel distribution! \n\n', idCellFull);
                else
                    fprintf('%s was not fitted with a kernel distribution! \n\n', idCellFull);
                end
                if isfield(params, 'mu1Model2')
                    ctCellModel2 = ctCellModel2 + 1;
                    sliceNumModel2(ctCellModel2, 1) = iSlice;
                    mu1Model2AllCells(ctCellModel2, 1) = params.mu1Model2;
                    mu2Model2AllCells(ctCellModel2, 1) = params.mu2Model2;
                    spacingModel2AllCells(ctCellModel2, 1) = params.spacingModel2;
                    threshold1Model2AllCells(ctCellModel2, 1) = params.threshold1Model2;
                    threshold2Model2AllCells(ctCellModel2, 1) = params.threshold2Model2;
                    void1Model2AllCells(ctCellModel2, 1) = params.void1Model2;
                    void2Model2AllCells(ctCellModel2, 1) = params.void2Model2;
                    if plotOverlapFlag || plotHeatMapFlag
                        pdfModel2AllCells{iSlice} = [pdfModel2AllCells{iSlice}; ...
                                                    params.pdfModel2(xVec)];
                        comp2Model2AllCells{iSlice} = [comp2Model2AllCells{iSlice}; ...
                                                    params.comp2Model2(xVec)];
                    end
                    fprintf('Done fitting %s with two Gaussians! \n\n', idCellFull);
                else
                    fprintf('%s cannot be fitted with two Gaussians! \n\n', ...
                            idCellFull);
                end

                if isfield(params, 'mu1Model3')
                    ctCellModel3 = ctCellModel3 + 1;
                    sliceNumModel3(ctCellModel3, 1) = iSlice;
                    mu1Model3AllCells(ctCellModel3, 1) = params.mu1Model3;
                    mu2Model3AllCells(ctCellModel3, 1) = params.mu2Model3;
                    spacingModel3AllCells(ctCellModel3, 1) = params.spacingModel3;
                    threshold1Model3AllCells(ctCellModel3, 1) = params.threshold1Model3;
                    threshold2Model3AllCells(ctCellModel3, 1) = params.threshold2Model3;
                    void1Model3AllCells(ctCellModel3, 1) = params.void1Model3;
                    void2Model3AllCells(ctCellModel3, 1) = params.void2Model3;
                    if plotOverlapFlag || plotHeatMapFlag
                        pdfModel3AllCells{iSlice} = [pdfModel3AllCells{iSlice}; ...
                                                    params.pdfModel3(xVec)];
                        comp2Model3AllCells{iSlice} = [comp2Model3AllCells{iSlice}; ...
                                                    params.comp2Model3(xVec)];
                    end
                    fprintf('Done fitting %s with Gaussian-Exp-Exponential! \n\n', ...
                            idCellFull);
                else
                    fprintf('%s cannot be fitted with Gaussian-Exp-Exponential! \n\n', ...
                            idCellFull);
                end
            else
                fprintf('No logIEIs for %s!\n\n', idCellFull);
            end
        end

        % Plot pdfs of all cells together
        if plotOverlapFlag
            if ~isempty(pdfModel0AllCells)
                h = figure('Visible', 'off');
                hold on;
                plot(xVec, pdfModel0AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Kernel distributions for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(outFolderParent, idExpt, ...
                        [idExpt, '_', idSlice, '_pdfs_Kernel']), 'png');
                close(h);
            end

            if ~isempty(pdfModel2AllCells)
                h = figure('Visible', 'off');
                hold on;
                plot(xVec, pdfModel2AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Two-Gaussian fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(outFolderParent, idExpt, ...
                        [idExpt, '_', idSlice, '_pdfs_GaussOnly']), 'png');
                close(h);

                h = figure('Visible', 'off');
                hold on;
                plot(xVec, comp2Model2AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['2nd Gaussian fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(outFolderParent, idExpt, ...
                        [idExpt, '_', idSlice, '_comp2s_GaussOnly']), 'png');
                close(h);
            end

            if ~isempty(pdfModel3AllCells)
                h = figure('Visible', 'off');
                hold on;
                plot(xVec, pdfModel3AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Gaussian-Exp-Exponential fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(outFolderParent, idExpt, ...
                        [idExpt, '_', idSlice, '_pdfs_GaussExp']), 'png');
                close(h);

                h = figure('Visible', 'off');
                hold on;
                plot(xVec, comp2Model3AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Exp-Exponential fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(outFolderParent, idExpt, ...
                        [idExpt, '_', idSlice, '_comp2s_GaussExp']), 'png');
                close(h);
            end
        end
    end
    % Extract information from fits
    commands = generate_stats_commands('Cell', 'Model0');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Cell', 'Model2');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Cell', 'Model3');
    cellfun(@eval, commands);

    for iSlice = 1:nSlices      % for each slice
        % Extract identifier for slice
        idSlice = idAllSlices{iSlice};
        fprintf('Slice %s ... \n', idSlice);

        % Extract IEIs for each cell of this slice and remove pooled logIEIs
        thisSlice = thisExpt.(idSlice).perCell;             % a cell array
        nCells = numel(thisSlice);
        
        if nCells > 0 && plotHeatMapFlag
            commands = generate_heatmap_commands('Cell', 'Model0', outFolderParent);
            cellfun(@eval, commands);
            commands = generate_heatmap_commands('Cell', 'Model2', outFolderParent);
            cellfun(@eval, commands);
            commands = generate_heatmap_commands('Cell', 'Model3', outFolderParent);
            cellfun(@eval, commands);
        end
    end

    % Plot histograms of statistics of all cells grouped by slice
    if plotStatHistFlag
        figname = fullfile(outFolderParent, ['mu1_', idExpt, '_Kernel']);
        plot_grouped_histogram(figname, exp(mu1Model0AllCells), ...
                    sliceNumModel0, idAllSlices, ...
                    'Intra-burst means', 's', ...
                    ['Intra-burst means of ', idExpt, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['mu1_', idExpt, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(mu1Model2AllCells), ...
                    sliceNumModel2, idAllSlices, ...
                    'Intra-burst means', 's', ...
                    ['Intra-burst means of ', idExpt, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['mu1_', idExpt, '_GaussExp']);
        plot_grouped_histogram(figname, exp(mu1Model3AllCells), ...
                    sliceNumModel3, idAllSlices, ...
                    'Intra-burst means', 's', ...
                    ['Intra-burst means of ', idExpt, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['mu2_', idExpt, '_Kernel']);
        plot_grouped_histogram(figname, exp(mu2Model0AllCells), ...
                    sliceNumModel0, idAllSlices, ...
                    'Inter-burst means', 's', ...
                    ['Inter-burst means of ', idExpt, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['mu2_', idExpt, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(mu2Model2AllCells), ...
                    sliceNumModel2, idAllSlices, ...
                    'Inter-burst means', 's', ...
                    ['Inter-burst means of ', idExpt, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['mu2_', idExpt, '_GaussExp']);
        plot_grouped_histogram(figname, mu2Model3AllCells, ...
                    sliceNumModel3, idAllSlices, ...
                    'Inter-burst means', 's', ...
                    ['Inter-burst means of ', idExpt, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['spacing_', idExpt, '_Kernel']);
        plot_grouped_histogram(figname, spacingModel0AllCells, ...
                    sliceNumModel0, idAllSlices, ...
                    'Spacing parameters', 'log s', ...
                    ['Spacing parameters of ', idExpt, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['spacing_', idExpt, '_GaussOnly']);
        plot_grouped_histogram(figname, spacingModel2AllCells, ...
                    sliceNumModel2, idAllSlices, ...
                    'Spacing parameters', 'log s', ...
                    ['Spacing parameters of ', idExpt, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['spacing_', idExpt, '_GaussExp']);
        plot_grouped_histogram(figname, spacingModel3AllCells, ...
                    sliceNumModel3, idAllSlices, ...
                    'Spacing parameters', 'log s', ...
                    ['Spacing parameters of ', idExpt, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['threshold_', idExpt, '_Kernel']);
        plot_grouped_histogram(figname, exp(threshold1Model0AllCells), ...
                    sliceNumModel0, idAllSlices, ...
                    'Thresholds', 's', ...
                    ['Thresholds of ', idExpt, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['threshold1_', idExpt, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(threshold1Model2AllCells), ...
                    sliceNumModel2, idAllSlices, ...
                    'Threshold #1s', 's', ...
                    ['Threshold #1s of ', idExpt, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['threshold1_', idExpt, '_GaussExp']);
        plot_grouped_histogram(figname, exp(threshold1Model3AllCells), ...
                    sliceNumModel3, idAllSlices, ...
                    'Threshold #1s', 's', ...
                    ['Threshold #1s of ', idExpt, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['threshold2_', idExpt, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(threshold2Model2AllCells), ...
                    sliceNumModel2, idAllSlices, ...
                    'Threshold #2s', 's', ...
                    ['Threshold #2s of ', idExpt, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['threshold2_', idExpt, '_GaussExp']);
        plot_grouped_histogram(figname, exp(threshold2Model3AllCells), ...
                    sliceNumModel3, idAllSlices, ...
                    'Threshold #2s', 's', ...
                    ['Threshold #2s of ', idExpt, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(outFolderParent, ['void_', idExpt, '_Kernel']);
        plot_grouped_histogram(figname, void1Model0AllCells, ...
                    sliceNumModel0, idAllSlices, ...
                    'Void parameters', '', ...
                    ['Void parameters of ', idExpt, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(outFolderParent, ['void1_', idExpt, '_GaussOnly']);
        plot_grouped_histogram(figname, void1Model2AllCells, ...
                    sliceNumModel2, idAllSlices, ...
                    'Void parameter #1s', '', ...
                    ['Void parameter #1s of ', idExpt, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(outFolderParent, ['void1_', idExpt, '_GaussExp']);
        plot_grouped_histogram(figname, void1Model3AllCells, ...
                    sliceNumModel3, idAllSlices, ...
                    'Void parameter #1s', '', ...
                    ['Void parameter #1s of ', idExpt, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(outFolderParent, ['void2_', idExpt, '_GaussOnly']);
        plot_grouped_histogram(figname, void2Model2AllCells, ...
                    sliceNumModel2, idAllSlices, ...
                    'Void parameter #2s', '', ...
                    ['Void parameter #2s of ', idExpt, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(outFolderParent, ['void2_', idExpt, '_GaussExp']);
        plot_grouped_histogram(figname, void2Model3AllCells, ...
                    sliceNumModel3, idAllSlices, ...
                    'Void parameter #2s', '', ...
                    ['Void parameter #2s of ', idExpt, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
    end

    % Plot scatter plots of statistics of all cells grouped by slice
    if plotScatterFlag
        commands = generate_scatter_commands(idExpt, 'Model0', outFolderParent);
        cellfun(@eval, commands);
        commands = generate_scatter_commands(idExpt, 'Model2', outFolderParent);
        cellfun(@eval, commands);
        commands = generate_scatter_commands(idExpt, 'Model3', outFolderParent);
        cellfun(@eval, commands);
    end
end

% Plot bar graphs comparing means
if plotBarFlag
    h = figure(1017);
    clf(h);
    h = bar_w_CI(h, muPerCellModel0, muPerCellModel0Low, ...
                    muPerCellModel0High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using kernel distributions for each cell');
    saveas(h, fullfile(outFolderParent, 'Kernel_means_AvgCells'), 'png');
    
    h = figure(1018);
    clf(h);
    h = bar_w_CI(h, muPerCellModel2, muPerCellModel2Low, ...
                    muPerCellModel2High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using two-Gaussian fits for each cell');
    saveas(h, fullfile(outFolderParent, 'GaussOnly_means_AvgCells'), 'png');

    h = figure(1019);
    clf(h);
    h = bar_w_CI(h, muPerCellModel3, muPerCellModel3Low, ...
                    muPerCellModel3High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using Gaussian-Exp-Exponential fits for each cell');
    saveas(h, fullfile(outFolderParent, 'GaussExp_means_AvgCells'), 'png');

    h = figure(1020);
    clf(h);
    h = bar_w_CI(h, [muPerCellModel0(1, :); muPerCellModel2(1, :); ...
                        muPerCellModel3(1, :)], ...
                    [muPerCellModel0Low(1, :); muPerCellModel2Low(1, :); ...
                        muPerCellModel3Low(1, :)], ...
                    [muPerCellModel0High(1, :); muPerCellModel2High(1, :); ...
                        muPerCellModel3High(1, :)], ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Mean inter-spike interval (s)');
    title('Mean intraburst ISIs averaged across cells');
    saveas(h, fullfile(outFolderParent, 'Mu1_AvgCells'), 'png');

    h = figure(1021);
    clf(h);
    h = bar_w_CI(h, spacingPerCell, spacingPerCellLow, spacingPerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Time scale spacing parameter');
    title('Spacing parameters averaged across cells');
    saveas(h, fullfile(outFolderParent, 'Spacing_AvgCells'), 'png');

    h = figure(1022);
    clf(h);
    h = bar_w_CI(h, threshold1PerCell, threshold1PerCellLow, threshold1PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #1 (s)');
    title('Burst detection threshold #1 averaged across cells');
    saveas(h, fullfile(outFolderParent, 'Threshold1_AvgCells'), 'png');

    h = figure(1023);
    clf(h);
    h = bar_w_CI(h, threshold2PerCell, threshold2PerCellLow, threshold2PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #2 (s)');
    title('Burst detection threshold #2 averaged across cells');
    saveas(h, fullfile(outFolderParent, 'Threshold2_AvgCells'), 'png');

    h = figure(1024);
    clf(h);
    h = bar_w_CI(h, void1PerCell, void1PerCellLow, void1PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #1');
    title('Void parameter #1s averaged across cells');
    saveas(h, fullfile(outFolderParent, 'Void1_AvgCells'), 'png');

    h = figure(1025);
    clf(h);
    h = bar_w_CI(h, void2PerCell, void2PerCellLow, void2PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allExpts, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #2');
    title('Void parameter #2s averaged across cells');
    saveas(h, fullfile(outFolderParent, 'Void2_AvgCells'), 'png');
end

%% Old method
%% Store cutoff
% data.events.cutoff.cutoff_seconds = cutOff ;

%% per experimental group
perExFeelds = fieldnames(data.perExpt.periods) ;
for i = 1:size(perExFeelds,1)
    allPeriods = sort(data.perExpt.periods.(perExFeelds{i}).allCatted) ;
  
    allPeriods(allPeriods==0) = [] ;
    
%     period20IDX = find(allPeriods>=cutOff,1) ;
%     allPeriodsLimits = allPeriods(1:period20IDX) ;

%     data.events.cutoff.periods.perExpt.(perExFeelds{i}) = allPeriods ;
    
    data.stats.(perExFeelds{i}).allPeriods.mean = mean(allPeriods) ;
    data.stats.(perExFeelds{i}).allPeriods.stdev = std(allPeriods) ;
    
%     data.stats.(perExFeelds{i}).cutoffPeriods.mean = mean(allPeriodsLimits) ;
%     data.stats.(perExFeelds{i}).cutoffPeriods.stdev = std(allPeriodsLimits) ;

    clear allPeriods allPeriodsConstrain
end

%% now for all combined experiments
allPeriods = sort(data.perEverything.periods) ;
allPeriods(allPeriods==0) = [] ;
% period20IDX = find(allPeriods>=cutOff,1) ;
allPeriodsLimits = allPeriods ;
% data.events.cutoff.periods.perEverything = allPeriodsLimits ;

data.stats.perEverything.allPeriods.periods = allPeriods ;
% data.stats.perEverything.cutoffPeriods.periods = allPeriodsLimits ;
data.stats.perEverything.allPeriods.mean = mean(allPeriods) ;
data.stats.perEverything.allPeriods.stdev = std(allPeriods) ;
% data.stats.perEverything.cutoffPeriods.mean = mean(allPeriodsLimits) ;
% data.stats.perEverything.cutoffPeriods.stdev = std(allPeriodsLimits) ;
% data.thresholds.allExperiments.threshold = data.thresholds.allExperiments.mean + 2*data.thresholds.allExperiments.stdev ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function commands = generate_stats_commands(perWhat, idModel)

switch idModel
case 'Model0'
    idMethod = 'Kernel';
    iModel = 1;
case 'Model2'
    idMethod = 'GaussOnly';
    iModel = 2;
case 'Model3'
    idMethod = 'GaussExpExp';
    iModel = 3;
otherwise
    error('ID of model not recognized!');
end

commands{1} = sprintf([...
    'meanMu1', idModel, ' = nanmean(mu1', idModel, 'All', perWhat, 's);\n', ...
    'meanMu2', idModel, ' = nanmean(mu2', idModel, 'All', perWhat, 's);\n', ...
    'stdErrMu1', idModel, ' = nanstderr(mu1', idModel, 'All', perWhat, 's);\n', ...
    'stdErrMu2', idModel, ' = nanstderr(mu2', idModel, 'All', perWhat, 's);\n', ...
    'muPer', perWhat, idModel, '(1, iExpt) = exp(meanMu1', idModel, ');\n', ...
    'muPer', perWhat, idModel, 'Low(1, iExpt) = exp(meanMu1', idModel, ' - 1.96*stdErrMu1', idModel, ');\n', ...
    'muPer', perWhat, idModel, 'High(1, iExpt) = exp(meanMu1', idModel, ' + 1.96*stdErrMu1', idModel, ');\n', ...
    ]);
switch idModel
case {'Model0', 'Model2'}
    commands{2} = sprintf([...
        'muPer', perWhat, idModel, '(2, iExpt) =  ', ...
            'exp(meanMu2', idModel, ');\n', ...
        'muPer', perWhat, idModel, 'Low(2, iExpt) =  ', ...
            'exp(meanMu2', idModel, ' - 1.96*stdErrMu2', idModel, ');\n', ...
        'muPer', perWhat, idModel, 'High(2, iExpt) =  ', ...
            'exp(meanMu2', idModel, ' + 1.96*stdErrMu2', idModel, ');\n', ...
        ]);
case 'Model3'
    commands{2} = sprintf([...
        'muPer', perWhat, idModel, '(2, iExpt) =  ', ...
            'meanMu2', idModel, ';\n', ...
        'muPer', perWhat, idModel, 'Low(2, iExpt) =  ', ...
            'meanMu2', idModel, ' - 1.96*stdErrMu2', idModel, ';\n', ...
        'muPer', perWhat, idModel, 'High(2, iExpt) =  ', ...
            'meanMu2', idModel, ' + 1.96*stdErrMu2', idModel, ';\n', ...
        ]);
otherwise
    error('ID of model not recognized!');
end
commands{3} = sprintf([...
    'meanSpacing', idModel, ' = nanmean(spacing', idModel, 'All', perWhat, 's);\n', ...
    'stdErrSpacing', idModel, ' = nanstderr(spacing', idModel, 'All', perWhat, 's);\n', ...
    'spacingPer', perWhat, '(', num2str(iModel), ', iExpt) = ', ...
		'meanSpacing', idModel, ';\n', ...
    'spacingPer', perWhat, 'Low(', num2str(iModel), ', iExpt) = ', ...
		'meanSpacing', idModel, ' - 1.96*stdErrSpacing', idModel, ';\n', ...
    'spacingPer', perWhat, 'High(', num2str(iModel), ', iExpt) = ', ...
		'meanSpacing', idModel, ' + 1.96*stdErrSpacing', idModel, ';\n', ...
    'meanThreshold1', idModel, ' = nanmean(threshold1', idModel, 'All', perWhat, 's);\n', ...
    'stdErrThreshold1', idModel, ' = nanstderr(threshold1', idModel, 'All', perWhat, 's);\n', ...
    'threshold1Per', perWhat, '(', num2str(iModel), ', iExpt) = ', ...
		'exp(meanThreshold1', idModel, ');\n', ...
    'threshold1Per', perWhat, 'Low(', num2str(iModel), ', iExpt) = ', ...
		'exp(meanThreshold1', idModel, ' - 1.96*stdErrThreshold1', idModel, ');\n', ...
    'threshold1Per', perWhat, 'High(', num2str(iModel), ', iExpt) = ', ...
		'exp(meanThreshold1', idModel, ' + 1.96*stdErrThreshold1', idModel, ');\n', ...
    'meanThreshold2', idModel, ' = nanmean(threshold2', idModel, 'All', perWhat, 's);\n', ...
    'stdErrThreshold2', idModel, ' = nanstderr(threshold2', idModel, 'All', perWhat, 's);\n', ...
    'threshold2Per', perWhat, '(', num2str(iModel), ', iExpt) = ', ...
		'exp(meanThreshold2', idModel, ');\n', ...
    'threshold2Per', perWhat, 'Low(', num2str(iModel), ', iExpt) = ', ...
		'exp(meanThreshold2', idModel, ' - 1.96*stdErrThreshold2', idModel, ');\n', ...
    'threshold2Per', perWhat, 'High(', num2str(iModel), ', iExpt) = ', ...
		'exp(meanThreshold2', idModel, ' + 1.96*stdErrThreshold2', idModel, ');\n', ...
    'meanVoid1', idModel, ' = nanmean(void1', idModel, 'All', perWhat, 's);\n', ...
    'stdErrVoid1', idModel, ' = nanstderr(void1', idModel, 'All', perWhat, 's);\n', ...
    'void1Per', perWhat, '(', num2str(iModel), ', iExpt) = ', ...
		'meanVoid1', idModel, ';\n', ...
    'void1Per', perWhat, 'Low(', num2str(iModel), ', iExpt) = ', ...
		'meanVoid1', idModel, ' - 1.96*stdErrVoid1', idModel, ';\n', ...
    'void1Per', perWhat, 'High(', num2str(iModel), ', iExpt) = ', ...
		'meanVoid1', idModel, ' + 1.96*stdErrVoid1', idModel, ';\n', ...
    'meanVoid2', idModel, ' = nanmean(void2', idModel, 'All', perWhat, 's);\n', ...
    'stdErrVoid2', idModel, ' = nanstderr(void2', idModel, 'All', perWhat, 's);\n', ...
    'void2Per', perWhat, '(', num2str(iModel), ', iExpt) = ', ...
		'meanVoid2', idModel, ';\n', ...
    'void2Per', perWhat, 'Low(', num2str(iModel), ', iExpt) = ', ...
		'meanVoid2', idModel, ' - 1.96*stdErrVoid2', idModel, ';\n', ...
    'void2Per', perWhat, 'High(', num2str(iModel), ', iExpt) = ', ...
		'meanVoid2', idModel, ' + 1.96*stdErrVoid2', idModel, ';\n', ...
    ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function commands = generate_heatmap_commands(perWhat, idModel, outFolderParent);

switch idModel
case 'Model0'
    idMethod = 'Kernel';
    descriptionMethod = 'Kernel distributions';
case 'Model2'
    idMethod = 'GaussOnly';
    descriptionMethod = 'Two-Gaussian distributions';
case 'Model3'
    idMethod = 'GaussExpExp';
    descriptionMethod = 'Gaussian-ExpExponential distributions';
otherwise
    error('ID of model not recognized!');
end

switch perWhat
case 'Slice'
    str1 = '';
    str2 = '';
    str3 = '';
    str4 = 'idExpt';
case 'Cell'
    str1 = '{iSlice}';
    str2 = ['(sliceNum', idModel, ' == iSlice)'];
    str3 = 'idExpt, ''/'', ';
    str4 = 'idExpt, ''_'', idSlice';
otherwise
    error('perWhat not recognized!');
end

commands{1} = sprintf([...
    'if ~isempty(pdf', idModel, 'All', perWhat, 's', str1, ')\n', ...
        'n', perWhat, 'sUsed = size(pdf', idModel, 'All', perWhat, 's', str1, ', 1);\n', ...
        'h = figure(''Visible'', ''off'');\n', ...
        'cm = colormap(hot);\n', ...
        'nColors = size(cm, 1);\n', ...
        '[~, origInd] = sort(void2', idModel, 'All', perWhat, 's', str2, ', ''descend'');\n', ...
        'cData = color_index(pdf', idModel, 'All', perWhat, 's', str1, ', ', ...
            'union(linspace(0, cutPdfValue, floor(nColors/2) + 1), ', ...
            'linspace(cutPdfValue, maxPdfValue, nColors - floor(nColors/2) + 1)));\n', ...
        'I = imagesc([min(xVec), max(xVec)], ', ...
            '[1, n', perWhat, 'sUsed], cData(origInd, :));\n', ...
        'set(I, ''CDataMapping'', ''direct'');\n', ...
        'set(gca, ''YTick'', 1:n', perWhat, 'sUsed + 0.5);\n', ...
        'set(gca, ''YTickLabel'', idAll', perWhat, 's', str1, '(origInd));\n', ...
        'xlabel(xUnitHist);\n', ...
        'colorbar(''Ticks'', [1, floor(nColors/4), floor(nColors/2), ', ...
            'floor(nColors*3/4), nColors], ', ...
            '''TickLabels'', [0, cutPdfValue/2, cutPdfValue, ', ...
            '(cutPdfValue + maxPdfValue)/2, maxPdfValue]);\n', ...
        'title([''', descriptionMethod, ' for ', perWhat, 's of '', ', ...
            str4, '], ''Interpreter'', ''none'');\n', ...
        'saveas(h, fullfile(''', outFolderParent, ''', [', str3, '''heatmap_'', ', str4, ...
        ', ''_pdfs_', idMethod, ''']), ''png'');\n', ...
        'close(h);\n', ...
    'end', ...
    ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function commands = generate_scatter_commands(idExpt, idModel, outFolderParent)

switch idModel
case 'Model0'
    idMethod = 'Kernel';
    descriptionMethod = 'kernel distribution fits';
    str1 = 'exp(mu2';
    str2 = 'AllCells), ';
case 'Model2'
    idMethod = 'GaussOnly';
    descriptionMethod = 'Two-Gaussian fits';
    str1 = 'exp(mu2';
    str2 = 'AllCells), ';
case 'Model3'
    idMethod = 'GaussExpExp';
    descriptionMethod = 'Gaussian-ExpExponential fits';
    str1 = 'mu2';
    str2 = 'AllCells, ';    
otherwise
    error('ID of model not recognized!');
end

figname = fullfile(outFolderParent, ['scatter_mu1_mu2_', idExpt, '_', idMethod, '_unscaled']);
commands{1} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(mu1', idModel, 'AllCells), ', str1, idModel, str2, ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Intra-burst mean'', ''s'', ', ...
            '''Inter-burst mean'', ''s'', ', ...
            '''Means for ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false);'];

figname = fullfile(outFolderParent, ['scatter_mu1_mu2_', idExpt, '_', idMethod, '_log']);
commands{2} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(mu1', idModel, 'AllCells), ', str1, idModel, str2, ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Intra-burst mean'', ''s'', ', ...
            '''Inter-burst mean'', ''s'', ', ...
            '''Means for ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''log'', ', ...
            '''XLimits'', [0.5, 8], ''YLimits'', [1, 1000]);'];

figname = fullfile(outFolderParent, ['scatter_mu1_mu2_', idExpt, '_', idMethod, '_linear']);
commands{3} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(mu1', idModel, 'AllCells), ', str1, idModel, str2, ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Intra-burst mean'', ''s'', ', ...
            '''Inter-burst mean'', ''s'', ', ...
            '''Means for ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false, ', ...
            '''XLimits'', [0.5, 7.5], ''YLimits'', [0, 220]);'];

figname = fullfile(outFolderParent, ['scatter_threshold1_threshold2_', idExpt, '_', idMethod, '_unscaled']);
commands{4} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold1', idModel, 'AllCells), exp(threshold2', idModel, 'AllCells), ', ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Threshold (intersection)'', ''s'', ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Thresholds for ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false);'];

figname = fullfile(outFolderParent, ['scatter_threshold1_threshold2_', idExpt, '_', idMethod, '_log']);
commands{5} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold1', idModel, 'AllCells), exp(threshold2', idModel, 'AllCells), ', ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Threshold (intersection)'', ''s'', ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Thresholds for ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''log'', ', ...
            '''XLimits'', [2, 70], ''YLimits'', [2, 200]);'];

figname = fullfile(outFolderParent, ['scatter_threshold1_threshold2_', idExpt, '_', idMethod, '_linear']);
commands{6} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold1', idModel, 'AllCells), exp(threshold2', idModel, 'AllCells), ', ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Threshold (intersection)'', ''s'', ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Thresholds for ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false, ', ...
            '''XLimits'', [0, 30], ''YLimits'', [0, 40]);'];

figname = fullfile(outFolderParent, ['scatter_void1_void2_', idExpt, '_', idMethod]);
commands{7} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'void1', idModel, 'AllCells, void2', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Void Parameter (intersection)'', '''', ', ...
            '''Void Parameter (minimum)'', '''', ', ...
            '''Void Parameters for ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false, ', ...
            '''XLimits'', [0, 1], ''YLimits'', [0, 1]);'];

figname = fullfile(outFolderParent, ['scatter_threshold2_void2_', idExpt, '_', idMethod]);
commands{8} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold2', idModel, 'AllCells), void2', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Void Parameter (minimum)'', '''', ', ...
            '''Void Parameter vs. Threshold for cells in ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''linear'', ', ...
            '''XLimits'', [2, 200], ''YLimits'', [0, 1]);'];

figname = fullfile(outFolderParent, ['scatter_threshold2_spacing_', idExpt, '_', idMethod]);
commands{9} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold2', idModel, 'AllCells), spacing', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Spacing Parameter'', '''', ', ...
            '''Spacing Parameter vs. Threshold for cells in ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''linear'', ', ...
            '''XLimits'', [2, 200], ''YLimits'', [0, 5]);'];

figname = fullfile(outFolderParent, ['scatter_spacing_void2_', idExpt, '_', idMethod]);
commands{10} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'spacing', idModel, 'AllCells, void2', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', idAllSlices, ', ...
            '''Spacing Parameter'', '''', ', ...
            '''Void Parameter (minimum)'', '''', ', ...
            '''Burstiness of cells in ', idExpt, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XLimits'', [0, 5], ''YLimits'', [0, 1]);'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:


    % Store stats in data structure
    data.stats.GaussExpExp.(idExpt).perSlice.meanMu1 = meanMu1Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.meanMu2 = meanMu2Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.stdErrMu1 = stdErrMu1Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.stdErrMu2 = stdErrMu2Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.meanSpacing = meanSpacingModel3;
    data.stats.GaussExpExp.(idExpt).perSlice.stdErrSpacing =stdErrSpacingModel3;
    data.stats.GaussExpExp.(idExpt).perSlice.meanThreshold1 = meanThreshold1Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.stdErrThreshold1 =stdErrThreshold1Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.meanThreshold2 = meanThreshold2Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.stdErrThreshold2 =stdErrThreshold2Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.meanVoid1 = meanVoid1Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.stdErrVoid1 =stdErrVoid1Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.meanVoid2 = meanVoid2Model3;
    data.stats.GaussExpExp.(idExpt).perSlice.stdErrVoid2 =stdErrVoid2Model3;

    % Store Model 3 stats in data structure
    data.stats.GaussExpExp.(idExpt).perCell.meanMu1 = meanMu1Model3;
    data.stats.GaussExpExp.(idExpt).perCell.meanMu2 = meanMu2Model3;
    data.stats.GaussExpExp.(idExpt).perCell.stdErrMu1 = stdErrMu1Model3;
    data.stats.GaussExpExp.(idExpt).perCell.stdErrMu2 =stdErrMu2Model3;
    data.stats.GaussExpExp.(idExpt).perCell.meanSpacing = meanSpacingModel3;
    data.stats.GaussExpExp.(idExpt).perCell.stdErrSpacing =stdErrSpacingModel3;
    data.stats.GaussExpExp.(idExpt).perCell.meanThreshold1 = meanThreshold1Model3;
    data.stats.GaussExpExp.(idExpt).perCell.stdErrThreshold1 =stdErrThreshold1Model3;
    data.stats.GaussExpExp.(idExpt).perCell.meanThreshold2 = meanThreshold2Model3;
    data.stats.GaussExpExp.(idExpt).perCell.stdErrThreshold2 =stdErrThreshold2Model3;
    data.stats.GaussExpExp.(idExpt).perCell.meanVoid1 = meanVoid1Model3;
    data.stats.GaussExpExp.(idExpt).perCell.stdErrVoid1 =stdErrVoid1Model3;
    data.stats.GaussExpExp.(idExpt).perCell.meanVoid2 = meanVoid2Model3;
    data.stats.GaussExpExp.(idExpt).perCell.stdErrVoid2 =stdErrVoid2Model3;

                                          if plotting

paramsByExpt = cell(nExpts, 1);
        % Store parameter fits in cell array
        paramsByExpt{iExpt} = params;
% Store parameter fits in data structure
for iExpt = 1:nExpts            % for each experiment
    fits.(idExpt) = paramsByExpt{idExpt};
end

paramsBySlice = cell(nSlices, 1);
        % Store parameter fits in cell array
        paramsBySlice{idSlice} = params;
% Store parameter fits in data structure
for iSlice = 1:nSlices
    % Store parameter fits in data structure
    fits.(idExpt).(idSlice) = paramsBySlice{idSlice};
end

paramsByCell = cell(nCells, 1);
        % Store parameter fits in cell array
        paramsByCell{iCell} = params;
% Store parameter fits in data structure
for iCell = 1:nCells
    fits.(idExpt).(idSlice).(idCell) = paramsByCell{iCell};
end

%}
