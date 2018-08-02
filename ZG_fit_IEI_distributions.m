function fitsGrouped = ZG_fit_IEI_distributions (ieisGrouped, varargin)
%% Fit inter-event-interval distributions and log distributions
% Usage: fitsGrouped = ZG_fit_IEI_distributions (ieisGrouped, varargin)
% Outputs:
%       TODO
% Arguments:
%       ieisGrouped - inter-event intervals grouped by experimental group,
%                       slice and cell
%                   must be a nonempty structure
%       varargin    - 'OutFolderParent': the output directory
%                   must be a valid directory
%                   default == pwd
%                   - 'FigsFolder': the figure output directory
%                   must be a valid directory
%                   default == fullfile(outFolderParent, 'FitFigures');
%                   - 'OutFileBase': the output file base name
%                   must be a string scalar or a character vector
%                   default == 'fitsGrouped'
%
% Requires:
%       /home/Matlab/Adams_Functions/fit_IEI.m
%       /home/Matlab/Adams_Functions/fit_logIEI.m
%       /home/Matlab/Adams_Functions/nanstderr.m
%       /home/Matlab/Adams_Functions/color_index.m
%       /home/Matlab/Adams_Functions/bar_w_CI.m
%       /home/Matlab/Adams_Functions/plot_grouped_histogram.m
%       /home/Matlab/Adams_Functions/plot_grouped_scatter.m
%
% Used by:
%       /home/Matlab/Adams_Functions/ZG_compute_IEI_thresholds.m
%
% File History:
% 2018-07-25 Adapted code from zgPeriodStats.m, 
%               which was derived from paula_iei4.m
% 2018-07-30 Now fits both IEIs and logIEIs
% TODO: Test and debug all plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
xUnitHist = 'log s';
plotFitsFlag = false;       % whether to plot individual fits
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

%% Set the order of groups
% allGroups = {'xHC067_3nMAng_Control', ...
%             'x10HC067_3nMAng', ...
%             'x50pMAng', 'x300pMAng', ...
%             'x3nMAng', 'x1uMang'};

%% Default values for optional arguments
outFolderParentDefault = pwd;
figsFolderDefault = '';
outFileBaseDefault = 'fitsGrouped';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'ieisGrouped', ...              % a structure
    @(x) validateattributes(x, {'struct'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutFolderParent', outFolderParentDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'FigsFolder', figsFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'OutFileBase', outFileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, ieisGrouped, varargin{:});
outFolderParent = iP.Results.OutFolderParent;
figsFolder = iP.Results.FigsFolder;
outFileBase = iP.Results.OutFileBase;

% Set defaults for dependent arguments
if isempty(figsFolder)
    figsFolder = fullfile(outFolderParent, 'FitFigures');
end

%% Create outFolderParent and figsFolder if they don't exist
if exist(outFolderParent, 'dir') ~= 7
    mkdir(outFolderParent);
end
if exist(figsFolder, 'dir') ~= 7
    mkdir(figsFolder);
end

%% Fit with distributions and compute thresholds for each experimental group
% Determine the experimental groups
allGroups = fieldnames(ieisGrouped);
nGroups = numel(allGroups);

% Fit and plot pooled IEIs and logIEIs for each experimental group
minIEIs = zeros(nGroups, 1);
maxIEIs = zeros(nGroups, 1);
minLogIEIs = zeros(nGroups, 1);
maxLogIEIs = zeros(nGroups, 1);
muModel0 = zeros(2, nGroups);
muModel2 = zeros(2, nGroups);
muModel2Low = zeros(2, nGroups);
muModel2High = zeros(2, nGroups);
muModel3 = zeros(2, nGroups);
muModel3Low = zeros(2, nGroups);
muModel3High = zeros(2, nGroups);
muModel4 = zeros(2, nGroups);
muModel4Low = zeros(2, nGroups);
muModel4High = zeros(2, nGroups);
muLogModel0 = zeros(2, nGroups);
muLogModel2 = zeros(2, nGroups);
muLogModel2Low = zeros(2, nGroups);
muLogModel2High = zeros(2, nGroups);
muLogModel3 = zeros(2, nGroups);
muLogModel3Low = zeros(2, nGroups);
muLogModel3High = zeros(2, nGroups);
for iGroup = 1:nGroups            % for each experimental group
    % Extract identifier for experimental group
    idGroup = allGroups{iGroup};
    fprintf('Processing group %s ... \n', idGroup);

    % Extract all IEIs pooled from different slices of this experimental group
    IEIsThisGroup = ieisGrouped.(idGroup).interEventIntervals;

    % Convert to log(IEI)s
    logIEIsThisGroup = log(IEIsThisGroup);

    % Determine maximum and minimum IEIs and logIEIs
    if ~isempty(IEIsThisGroup)
        minIEIs(iGroup) = min(IEIsThisGroup);
        maxIEIs(iGroup) = max(IEIsThisGroup);
    else
        minIEIs(iGroup) = NaN;
        maxIEIs(iGroup) = NaN; 
    end
    minLogIEIs(iGroup) = log(minIEIs(iGroup));
    maxLogIEIs(iGroup) = log(maxIEIs(iGroup)); 
 
    % Fit IEIs to distributions and plot on top of histograms
    if ~isempty(IEIsThisGroup)
        [~, ~, ~, ~, ~, IEIparams] = ...
            fit_IEI(IEIsThisGroup, idGroup, ...
                    'XUnit', xUnitHist, ...
                    'TruncateFlag', truncateDataFlag, ...
                    'PlotFlag', plotFitsFlag, ...
                    'OutFolder', figsFolder);
        fprintf('Done fitting IEIs of experimental group %s!\n', idGroup);

        % Store parameter fits in data structure
        fitsGrouped.(idGroup).IEIparams = IEIparams;

        % Extract means from fits
        if isfield(IEIparams, 'mu1Model0')
            muModel0(1, iGroup) = IEIparams.mu1Model0;
            muModel0(2, iGroup) = IEIparams.mu2Model0;
        else
            muModel0(1, iGroup) = NaN;
            muModel0(2, iGroup) = NaN;
        end
        if isfield(IEIparams, 'mu1Model2')
            muModel2(1, iGroup) = IEIparams.mu1Model2;
            muModel2Low(1, iGroup) = IEIparams.mu1CIModel2(1);
            muModel2High(1, iGroup) = IEIparams.mu1CIModel2(2);
            muModel2(2, iGroup) = IEIparams.mu2Model2;
            muModel2Low(2, iGroup) = IEIparams.mu2CIModel2(1);
            muModel2High(2, iGroup) = IEIparams.mu2CIModel2(2);
        else
            muModel2(1, iGroup) = NaN;
            muModel2Low(1, iGroup) = NaN;
            muModel2High(1, iGroup) = NaN;
            muModel2(2, iGroup) = NaN;
            muModel2Low(2, iGroup) = NaN;
            muModel2High(2, iGroup) = NaN;
        end
        if isfield(IEIparams, 'mu1Model3')
            muModel3(1, iGroup) = IEIparams.mu1Model3;
            muModel3Low(1, iGroup) = IEIparams.mu1CIModel3(1);
            muModel3High(1, iGroup) = IEIparams.mu1CIModel3(2);
            muModel3(2, iGroup) = IEIparams.mu2Model3;
            muModel3Low(2, iGroup) = IEIparams.mu2CIModel3(1);
            muModel3High(2, iGroup) = IEIparams.mu2CIModel3(2);
        else
            muModel3(1, iGroup) = NaN;
            muModel3Low(1, iGroup) = NaN;
            muModel3High(1, iGroup) = NaN;
            muModel3(2, iGroup) = NaN;
            muModel3Low(2, iGroup) = NaN;
            muModel3High(2, iGroup) = NaN;
        end
        if isfield(IEIparams, 'mean1Model4')
            muModel4(1, iGroup) = IEIparams.mean1Model4;
            muModel4Low(1, iGroup) = IEIparams.mean1CIModel4(1);
            muModel4High(1, iGroup) = IEIparams.mean1CIModel4(2);
            muModel4(2, iGroup) = IEIparams.mean2Model4;
            muModel4Low(2, iGroup) = IEIparams.mean2CIModel4(1);
            muModel4High(2, iGroup) = IEIparams.mean2CIModel4(2);
        else
            muModel4(1, iGroup) = NaN;
            muModel4Low(1, iGroup) = NaN;
            muModel4High(1, iGroup) = NaN;
            muModel4(2, iGroup) = NaN;
            muModel4Low(2, iGroup) = NaN;
            muModel4High(2, iGroup) = NaN;
        end
    else
        fprintf('No IEIs for experimental group %s!\n\n', idGroup);
    end

    % Fit logIEIs to distributions and plot on top of histograms
    if ~isempty(logIEIsThisGroup)
        [~, ~, ~, ~, logIEIparams] = ...
            fit_logIEI(logIEIsThisGroup, idGroup, ...
                    'XUnit', xUnitHist, ...
                    'TruncateFlag', truncateDataFlag, ...
                    'PlotFlag', plotFitsFlag, ...
                    'OutFolder', figsFolder);
        fprintf('Done fitting log(IEI)s of experimental group %s!\n', idGroup);

        % Store parameter fits in data structure
        fitsGrouped.(idGroup).logIEIparams = logIEIparams;

        % Extract means from fits
        if isfield(logIEIparams, 'mu1Model0')
            muModel0(1, iGroup) = exp(logIEIparams.mu1Model0);
            muModel0(2, iGroup) = exp(logIEIparams.mu2Model0);
        else
            muModel0(1, iGroup) = NaN;
            muModel0(2, iGroup) = NaN;
        end
        if isfield(logIEIparams, 'mu1Model2')
            muModel2(1, iGroup) = exp(logIEIparams.mu1Model2);
            muModel2Low(1, iGroup) = exp(logIEIparams.mu1CIModel2(1));
            muModel2High(1, iGroup) = exp(logIEIparams.mu1CIModel2(2));
            muModel2(2, iGroup) = exp(logIEIparams.mu2Model2);
            muModel2Low(2, iGroup) = exp(logIEIparams.mu2CIModel2(1));
            muModel2High(2, iGroup) = exp(logIEIparams.mu2CIModel2(2));
        else
            muModel2(1, iGroup) = NaN;
            muModel2Low(1, iGroup) = NaN;
            muModel2High(1, iGroup) = NaN;
            muModel2(2, iGroup) = NaN;
            muModel2Low(2, iGroup) = NaN;
            muModel2High(2, iGroup) = NaN;
        end
        if isfield(logIEIparams, 'mu1Model3')
            muModel3(1, iGroup) = exp(logIEIparams.mu1Model3);
            muModel3Low(1, iGroup) = exp(logIEIparams.mu1CIModel3(1));
            muModel3High(1, iGroup) = exp(logIEIparams.mu1CIModel3(2));
            muModel3(2, iGroup) = logIEIparams.mu2Model3;
            muModel3Low(2, iGroup) = logIEIparams.mu2CIModel3(1);
            muModel3High(2, iGroup) = logIEIparams.mu2CIModel3(2);
        else
            muModel3(1, iGroup) = NaN;
            muModel3Low(1, iGroup) = NaN;
            muModel3High(1, iGroup) = NaN;
            muModel3(2, iGroup) = NaN;
            muModel3Low(2, iGroup) = NaN;
            muModel3High(2, iGroup) = NaN;
        end
    else
        fprintf('No logIEIs for experimental group %s!\n\n', idGroup);
    end
end

% Compute x limits for plots; each row is an experimental group
xlimits = [minLogIEIs, maxLogIEIs];

% Compute x vectors for plots; each row is an experimental group
xVecs = cell(nGroups, 1);
parfor iGroup = 1:nGroups
    xVecs{iGroup} = linspace(xlimits(iGroup, 1), xlimits(iGroup, 2), nPoints);
end

% Plot bar graphs comparing means
if plotBarFlag
    h = figure(1004);
    clf(h);
    bar(muLogModel0);
    set(gca, 'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using a kernel distribution for pooled logIEIs');
    saveas(h, fullfile(figsFolder, 'Kernel_means_Pooled'), 'png');

    h = figure(1005);
    clf(h);
    h = bar_w_CI(h, muLogModel2, muLogModel2Low, muLogModel2High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using two-Gaussian fits for pooled logIEIs');
    saveas(h, fullfile(figsFolder, 'GaussOnly_means_Pooled'), 'png');

    h = figure(1006);
    clf(h);
    h = bar_w_CI(h, muLogModel3, muLogModel3Low, muLogModel3High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using Gaussian-Exp-Exponential fits for pooled logIEIs');
    saveas(h, fullfile(figsFolder, 'GaussExpExp_means_Pooled'), 'png');

    h = figure(1007);
    clf(h);
    h = bar_w_CI(h, [muLogModel2(1, :); muLogModel3(1, :)], ...
                    [muLogModel2Low(1, :); muLogModel3Low(1, :)], ...
                    [muLogModel2High(1, :); muLogModel3High(1, :)], ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Mean inter-spike interval (s)');
    title('Mean intraburst ISIs for pooled logIEIs');
    saveas(h, fullfile(figsFolder, 'Mu1_Pooled'), 'png');
end

%% Fit with distributions and compute thresholds for each slice

% Fit and plot logIEIs from each slice for each experimental group
muPerSliceLogModel0 = zeros(2, nGroups);    % means averaged from each slice
muPerSliceLogModel0Low = zeros(2, nGroups);
muPerSliceLogModel0High = zeros(2, nGroups);
muPerSliceLogModel2 = zeros(2, nGroups);    % means averaged from each slice
muPerSliceLogModel2Low = zeros(2, nGroups);
muPerSliceLogModel2High = zeros(2, nGroups);
muPerSliceLogModel3 = zeros(2, nGroups);    % means averaged from each slice
muPerSliceLogModel3Low = zeros(2, nGroups);
muPerSliceLogModel3High = zeros(2, nGroups);
spacingPerSlice = zeros(2, nGroups);     % spacing parameters averaged from each slice
spacingPerSliceLow = zeros(2, nGroups);
spacingPerSliceHigh = zeros(2, nGroups);
threshold1PerSlice = zeros(2, nGroups);  % threshold #1s averaged from each slice
threshold1PerSliceLow = zeros(2, nGroups);
threshold1PerSliceHigh = zeros(2, nGroups);
threshold2PerSlice = zeros(2, nGroups);  % threshold #2s averaged from each slice
threshold2PerSliceLow = zeros(2, nGroups);
threshold2PerSliceHigh = zeros(2, nGroups);
void1PerSlice = zeros(2, nGroups);       % void parameter #1s averaged from each slice
void1PerSliceLow = zeros(2, nGroups);
void1PerSliceHigh = zeros(2, nGroups);
void2PerSlice = zeros(2, nGroups);       % void parameter #2s averaged from each slice
void2PerSliceLow = zeros(2, nGroups);
void2PerSliceHigh = zeros(2, nGroups);
for iGroup = 1:nGroups            % for each experimental group
    % Extract identifier for experimental group
    idGroup = allGroups{iGroup};
    if plotOverlapFlag || plotHeatMapFlag
        xVec = xVecs{iGroup};
    end
    fprintf('Processing all slices of Group %s ... \n', idGroup);

    % Extract IEIs for this experimental group and remove pooled logIEIs
    thisGroup = rmfield(ieisGrouped.(idGroup), 'interEventIntervals');  
                                        % a structure

    % Get all slices
    allSliceLabels = fieldnames(thisGroup);
    nSlices = numel(allSliceLabels);
    ctSliceLogModel0 = 0;                  % counts useful slices for LogModel 0
    ctSliceLogModel2 = 0;                  % counts useful slices for LogModel 2
    ctSliceLogModel3 = 0;                  % counts useful slices for LogModel 3
    mu1LogModel0AllSlices = [];            % stores means of intra-burst intervals
    mu2LogModel0AllSlices = [];            % stores means of inter-burst intervals
    spacingLogModel0AllSlices = [];        % stores spacing parameters
    threshold1LogModel0AllSlices = [];     % stores burst detection threshold #1s
    threshold2LogModel0AllSlices = [];     % stores burst detection threshold #2s
    void1LogModel0AllSlices = [];          % stores void parameter #1s
    void2LogModel0AllSlices = [];          % stores void parameter #2s
    mu1LogModel2AllSlices = [];            % stores means of intra-burst intervals
    mu2LogModel2AllSlices = [];            % stores means of inter-burst intervals
    spacingLogModel2AllSlices = [];        % stores spacing parameters
    threshold1LogModel2AllSlices = [];     % stores burst detection threshold #1s
    threshold2LogModel2AllSlices = [];     % stores burst detection threshold #2s
    void1LogModel2AllSlices = [];          % stores void parameter #1s
    void2LogModel2AllSlices = [];          % stores void parameter #2s
    mu1LogModel3AllSlices = [];            % stores means of intra-burst intervals
    mu2LogModel3AllSlices = [];            % stores means of inter-burst intervals
    spacingLogModel3AllSlices = [];        % stores spacing parameters
    threshold1LogModel3AllSlices = [];     % stores burst detection threshold #1s
    threshold2LogModel3AllSlices = [];     % stores burst detection threshold #2s
    void1LogModel3AllSlices = [];          % stores void parameter #1s
    void2LogModel3AllSlices = [];          % stores void parameter #2s
    if plotOverlapFlag || plotHeatMapFlag
        pdfLogModel0AllSlices = [];        % stores pdfs for LogModel 0
        pdfLogModel2AllSlices = [];        % stores pdfs for LogModel 2
        comp2LogModel2AllSlices = [];      % stores component #2 pdfs for LogModel 2
        pdfLogModel3AllSlices = [];        % stores pdfs for LogModel 3
        comp2LogModel3AllSlices = [];      % stores component #2 pdfs for LogModel 3
    end

    % Construct outFolder for this experimental group
    outFolder = fullfile(figsFolder, idGroup);
    if exist(outFolder, 'dir') ~= 7 && ...
        (plotBarFlag  || plotFitsFlag || plotOverlapFlag || plotHeatMapFlag)
        mkdir(outFolder);
    end

    for iSlice = 1:nSlices      % for each slice
        % Extract identifier for slice
        idSlice = allSliceLabels{iSlice};
        idSliceFull = [idGroup, '_', idSlice];
        fprintf('Processing Slice %s ... \n', idSlice);

        % Extract IEIs pooled from all cells of this slice
        IEIsThisSlice = thisGroup.(idSlice).interEventIntervals;

        % Convert into log(IEI)s
        logIEIsThisSlice = log(IEIsThisSlice);

        % Fit IEIs to distributions and plot on top of histograms
        if ~isempty(IEIsThisSlice)
            [~, ~, ~, ~, ~, IEIparams] = ...
                fit_IEI(IEIsThisSlice, idSliceFull, ...
                        'XUnit', xUnitHist, ...
                        'TruncateFlag', truncateDataFlag, ...
                        'PlotFlag', plotFitsFlag, ...
                        'OutFolder', outFolder);

            % Store parameter fits in data structure
            fitsGrouped.(idGroup).(idSlice).IEIparams = IEIparams;
        end

        % Fit logIEIs to distributions and plot on top of histograms
        if ~isempty(logIEIsThisSlice)
            [~, ~, ~, ~, params] = ...
                fit_logIEI(logIEIsThisSlice, idSliceFull, ...
                        'XUnit', xUnitHist, ...
                        'TruncateFlag', truncateDataFlag, ...
                        'PlotFlag', plotFitsFlag, ...
                        'OutFolder', outFolder);

            % Store parameter fits in data structure
            fitsGrouped.(idGroup).(idSlice).logIEIparams = params;

            if isfield(params, 'mu2Model0')
                ctSliceLogModel0 = ctSliceLogModel0 + 1;
                mu1LogModel0AllSlices(ctSliceLogModel0, 1) = params.mu1Model0;
                mu2LogModel0AllSlices(ctSliceLogModel0, 1) = params.mu2Model0;
                spacingLogModel0AllSlices(ctSliceLogModel0, 1) = params.spacingModel0;
                threshold1LogModel0AllSlices(ctSliceLogModel0, 1) = params.thresholdModel0;
                threshold2LogModel0AllSlices(ctSliceLogModel0, 1) = params.thresholdModel0;
                void1LogModel0AllSlices(ctSliceLogModel0, 1) = params.voidModel0;
                void2LogModel0AllSlices(ctSliceLogModel0, 1) = params.voidModel0;
                if plotOverlapFlag || plotHeatMapFlag
                    pdfLogModel0AllSlices = [pdfLogModel0AllSlices; ...
                                            params.pdfModel0(xVec)];
                end
                fprintf('Done fitting %s with a kernel distribution! \n\n', idSliceFull);
            else
                fprintf('%s was not fitted with a kernel distribution! \n\n', idSliceFull);
            end
            if isfield(params, 'mu1Model2')
                ctSliceLogModel2 = ctSliceLogModel2 + 1;
                mu1LogModel2AllSlices(ctSliceLogModel2, 1) = params.mu1Model2;
                mu2LogModel2AllSlices(ctSliceLogModel2, 1) = params.mu2Model2;
                spacingLogModel2AllSlices(ctSliceLogModel2, 1) = params.spacingModel2;
                threshold1LogModel2AllSlices(ctSliceLogModel2, 1) = params.threshold1Model2;
                threshold2LogModel2AllSlices(ctSliceLogModel2, 1) = params.threshold2Model2;
                void1LogModel2AllSlices(ctSliceLogModel2, 1) = params.void1Model2;
                void2LogModel2AllSlices(ctSliceLogModel2, 1) = params.void2Model2;
                if plotOverlapFlag || plotHeatMapFlag
                    pdfLogModel2AllSlices = [pdfLogModel2AllSlices; ...
                                            params.pdfModel2(xVec)];
                    comp2LogModel2AllSlices = [comp2LogModel2AllSlices; ...
                                            params.comp2Model2(xVec)];
                end
                fprintf('Done fitting %s with two Gaussians! \n\n', idSliceFull);
            else
                fprintf('%s was not fitted with two Gaussians! \n\n', idSliceFull);
            end
            if isfield(params, 'mu1Model3')
                ctSliceLogModel3 = ctSliceLogModel3 + 1;
                mu1LogModel3AllSlices(ctSliceLogModel3, 1) = params.mu1Model3;
                mu2LogModel3AllSlices(ctSliceLogModel3, 1) = params.mu2Model3;
                spacingLogModel3AllSlices(ctSliceLogModel3, 1) = params.spacingModel3;
                threshold1LogModel3AllSlices(ctSliceLogModel3, 1) = params.threshold1Model3;
                threshold2LogModel3AllSlices(ctSliceLogModel3, 1) = params.threshold2Model3;
                void1LogModel3AllSlices(ctSliceLogModel3, 1) = params.void1Model3;
                void2LogModel3AllSlices(ctSliceLogModel3, 1) = params.void2Model3;
                if plotOverlapFlag || plotHeatMapFlag
                    pdfLogModel3AllSlices = [pdfLogModel3AllSlices; ...
                                            params.pdfModel3(xVec)];
                    comp2LogModel3AllSlices = [comp2LogModel3AllSlices; ...
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

    % Combine stats
    % commands = generate_stats_commands('Slice', 'Model0');
    % cellfun(@eval, commands);
    % commands = generate_stats_commands('Slice', 'Model2');
    % cellfun(@eval, commands);
    % commands = generate_stats_commands('Slice', 'Model3');
    % cellfun(@eval, commands);
    % commands = generate_stats_commands('Slice', 'Model4');
    % cellfun(@eval, commands);
    commands = generate_stats_commands('Slice', 'LogModel0');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Slice', 'LogModel2');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Slice', 'LogModel3');
    cellfun(@eval, commands);

    % Plot pdfs of all slices together
    if plotOverlapFlag
        if ~isempty(pdfLogModel0AllSlices)
            h = figure('Visible', 'off');
            hold on;
            plot(xVec, pdfLogModel0AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Kernel distributions for slices of ', idGroup], ...
                    'Interpreter', 'none');
            legend(allSliceLabels, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(figsFolder, [idGroup, '_pdfs_Kernel']), 'png');
            close(h);
        end

        if ~isempty(pdfLogModel2AllSlices)
            h = figure('Visible', 'off');
            hold on;
            plot(xVec, pdfLogModel2AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Two-Gaussian fits for slices of ', idGroup], ...
                    'Interpreter', 'none');
            legend(allSliceLabels, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(figsFolder, [idGroup, '_pdfs_GaussOnly']), 'png');
            close(h);

            h = figure('Visible', 'off');
            hold on;
            plot(xVec, comp2LogModel2AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['2nd Gaussian fits for slices of ', idGroup], ...
                    'Interpreter', 'none');
            legend(allSliceLabels, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(figsFolder, [idGroup, '_comp2s_GaussOnly']), 'png');
            close(h);
        end

        if ~isempty(pdfLogModel3AllSlices)
            h = figure('Visible', 'off');
            hold on;
            plot(xVec, pdfLogModel3AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Gaussian-Exp-Exponential fits for slices of ', idGroup], ...
                    'Interpreter', 'none');
            legend(allSliceLabels, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(figsFolder, [idGroup, '_pdfs_GaussExpExp']), 'png');
            close(h);

            h = figure('Visible', 'off');
            hold on;
            plot(xVec, comp2LogModel3AllSlices);
            xlabel(xUnitHist);
            ylabel('Probability density');
            title(['Exp-Exponential fits for slices of ', idGroup], ...
                    'Interpreter', 'none');
            legend(allSliceLabels, 'Interpreter', 'none', 'location', 'northwest');
            saveas(h, fullfile(figsFolder, [idGroup, '_comp2s_GaussExpExp']), 'png');
            close(h);
        end
    end
    if plotHeatMapFlag
        commands = generate_heatmap_commands('Slice', 'LogModel0', figsFolder);
        cellfun(@eval, commands);
        commands = generate_heatmap_commands('Slice', 'LogModel2', figsFolder);
        cellfun(@eval, commands);
        commands = generate_heatmap_commands('Slice', 'LogModel3', figsFolder);
        cellfun(@eval, commands);
    end
end

% Plot bar graph comparing means
if plotBarFlag
    h = figure(1007);
    clf(h);
    h = bar_w_CI(h, muPerSliceLogModel2, muPerSliceLogModel2Low, ...
                    muPerSliceLogModel2High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using two Gaussians for each slice');
    saveas(h, fullfile(figsFolder, 'GaussOnly_means_AvgSlices'), 'png');

    h = figure(1008);
    clf(h);
    h = bar_w_CI(h, muPerSliceLogModel3, muPerSliceLogModel3Low, ...
                    muPerSliceLogModel3High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using Gaussian-Exp-Exponential fits for each slice');
    saveas(h, fullfile(figsFolder, 'GaussExpExp_means_AvgSlices'), 'png');

    h = figure(1009);
    clf(h);
    h = bar_w_CI(h, [muPerSliceLogModel2(1, :); muPerSliceLogModel3(1, :)], ...
                    [muPerSliceLogModel2Low(1, :); muPerSliceLogModel3Low(1, :)], ...
                    [muPerSliceLogModel2High(1, :); muPerSliceLogModel3High(1, :)], ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Mean inter-spike interval (s)');
    title('Mean intraburst ISIs averaged across slices');
    saveas(h, fullfile(figsFolder, 'mu1_AvgSlices'), 'png');

    h = figure(1010);
    clf(h);
    h = bar_w_CI(h, spacingPerSlice, spacingPerSliceLow, spacingPerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Time scale spacing parameter');
    title('Spacing parameters averaged across slices');
    saveas(h, fullfile(figsFolder, 'Spacing_AvgSlices'), 'png');

    h = figure(1011);
    clf(h);
    h = bar_w_CI(h, threshold1PerSlice, threshold1PerSliceLow, threshold1PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #1 (s)');
    title('Burst detection threshold #1 averaged across slices');
    saveas(h, fullfile(figsFolder, 'Threshold1_AvgSlices'), 'png');

    h = figure(1012);
    clf(h);
    h = bar_w_CI(h, threshold2PerSlice, threshold2PerSliceLow, threshold2PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #2 (s)');
    title('Burst detection threshold #2 averaged across slices');
    saveas(h, fullfile(figsFolder, 'Threshold2_AvgSlices'), 'png');

    h = figure(1013);
    clf(h);
    h = bar_w_CI(h, void1PerSlice, void1PerSliceLow, void1PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #1');
    title('Void parameter #1s averaged across slices');
    saveas(h, fullfile(figsFolder, 'Void1_AvgSlices'), 'png');

    h = figure(1014);
    clf(h);
    h = bar_w_CI(h, void2PerSlice, void2PerSliceLow, void2PerSliceHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #2');
    title('Void parameter #2s averaged across slices');
    saveas(h, fullfile(figsFolder, 'Void2_AvgSlices'), 'png');
end

%% Fit with distributions and compute thresholds for each cell

% Fit and plot logIEIs from each cell for each experimental group
muPerCellLogModel0 = zeros(2, nGroups);     % means averaged from each cell
muPerCellLogModel0Low = zeros(2, nGroups);
muPerCellLogModel0High = zeros(2, nGroups);
muPerCellLogModel2 = zeros(2, nGroups);     % means averaged from each cell
muPerCellLogModel2Low = zeros(2, nGroups);
muPerCellLogModel2High = zeros(2, nGroups);
muPerCellLogModel3 = zeros(2, nGroups);     % means averaged from each cell
muPerCellLogModel3Low = zeros(2, nGroups);
muPerCellLogModel3High = zeros(2, nGroups);
spacingPerCell = zeros(2, nGroups);      % spacing parameters averaged from each cell
spacingPerCellLow = zeros(2, nGroups);
spacingPerCellHigh = zeros(2, nGroups);
threshold1PerCell = zeros(2, nGroups);    % threshold #1s averaged from each cell
threshold1PerCellLow = zeros(2, nGroups);
threshold1PerCellHigh = zeros(2, nGroups);
void1PerCell = zeros(2, nGroups);         % void parameter #1s averaged from each cell
void1PerCellLow = zeros(2, nGroups);
void1PerCellHigh = zeros(2, nGroups);
for iGroup = 1:nGroups            % for each experimental group
    % Extract identifier for experimental group
    idGroup = allGroups{iGroup};
    if plotOverlapFlag || plotHeatMapFlag
        xVec = xVecs{iGroup};
    end
    fprintf('Processing all cells of Group %s ... \n', idGroup);

    % Extract IEIs for this experimental group and remove pooled logIEIs
    thisGroup = rmfield(ieisGrouped.(idGroup), 'interEventIntervals');
                                        % a structure

    % Get all slices
    allSliceLabels = fieldnames(thisGroup);
    nSlices = numel(allSliceLabels);
    ctCellLogModel0 = 0;               % counts useful cells
    sliceNumLogModel0 = [];            % stores slice number of all cells
    mu1LogModel0AllCells = [];         % stores means of intra-burst intervals
    mu2LogModel0AllCells = [];         % stores means of inter-burst intervals
    spacingLogModel0AllCells = [];     % stores spacing parameters
    threshold1LogModel0AllCells = [];  % stores burst detection threshold #1s
    threshold2LogModel0AllCells = [];  % stores burst detection threshold #2s
    void1LogModel0AllCells = [];       % stores void parameter #1s
    void2LogModel0AllCells = [];       % stores void parameter #2s
    ctCellLogModel2 = 0;               % counts useful cells
    sliceNumLogModel2 = [];            % stores slice number of all cells
    mu1LogModel2AllCells = [];         % stores means of intra-burst intervals
    mu2LogModel2AllCells = [];         % stores means of inter-burst intervals
    spacingLogModel2AllCells = [];     % stores spacing parameters
    threshold1LogModel2AllCells = [];  % stores burst detection threshold #1s
    threshold2LogModel2AllCells = [];  % stores burst detection threshold #2s
    void1LogModel2AllCells = [];       % stores void parameter #1s
    void2LogModel2AllCells = [];       % stores void parameter #2s
    ctCellLogModel3 = 0;               % counts useful cells
    sliceNumLogModel3 = [];            % stores slice number of all cells
    mu1LogModel3AllCells = [];         % stores means of intra-burst intervals
    mu2LogModel3AllCells = [];         % stores means of inter-burst intervals
    spacingLogModel3AllCells = [];     % stores spacing parameters
    threshold1LogModel3AllCells = [];  % stores burst detection threshold #1s
    threshold2LogModel3AllCells = [];  % stores burst detection threshold #2s
    void1LogModel3AllCells = [];       % stores void parameter #1s
    void2LogModel3AllCells = [];       % stores void parameter #2s
    if plotOverlapFlag || plotHeatMapFlag
        pdfLogModel0AllCells = cell(1, nSlices);     % stores pdfs for LogModel 0
        pdfLogModel2AllCells = cell(1, nSlices);     % stores pdfs for LogModel 2
        comp2LogModel2AllCells = cell(1, nSlices);   % stores component #2 pdfs for LogModel 2
        pdfLogModel3AllCells = cell(1, nSlices);     % stores pdfs for LogModel 3
        comp2LogModel3AllCells = cell(1, nSlices);   % stores component #2 pdfs for LogModel 3

        idAllCells = cell(nSlices); % stores IDs of all cells in the slice
    end
    for iSlice = 1:nSlices      % for each slice
        % Extract identifier for slice
        idSlice = allSliceLabels{iSlice};
        fprintf('Processing Slice %s ... \n', idSlice);

        % Construct outFolder for this slice
        outFolder = fullfile(figsFolder, idGroup, idSlice);
        if exist(outFolder, 'dir') ~= 7 && (plotBarFlag || plotFitsFlag)
            mkdir(outFolder);
        end

        % Extract IEIs for each cell of this slice and remove pooled logIEIs
        thisSlice = rmfield(thisGroup.(idSlice), 'interEventIntervals');
                                            % a structure

        % Get all cells
        allCellLabels = fieldnames(thisSlice);
        nCells = numel(allCellLabels);
        
        for iCell = 1:nCells    % for each cell
            % Construct identifier for this cell
            idCell = allCellLabels{iCell};
            idAllCells{iSlice}{iCell} = idCell;
            idCellFull = [idGroup, '_', idSlice, '_', idCell];
            fprintf('Processing Cell #%d ... \n', iCell);
            
            % Extract IEIs for this cell
            IEIsThisCell = thisSlice.(idCell).interEventIntervals;

            % Convert to log(IEI)s
            logIEIsThisCell = log(IEIsThisCell);

            % Fit logIEIs to distributions and plot on top of histograms
            if ~isempty(IEIsThisCell) && size(IEIsThisCell, 1) > 1
                [~, ~, ~, ~, ~, IEIparams] = ...
                    fit_IEI(IEIsThisCell, idCellFull, 'XUnit', xUnitHist, ...
                            'TruncateFlag', truncateDataFlag, ...
                            'PlotFlag', all([plotFitsFlag, plotPerCellFlag]), ...
                            'OutFolder', outFolder);

                % Store parameter fits in data structure
                fitsGrouped.(idGroup).(idSlice).(idCell).IEIparams = IEIparams;
            end

            % Fit logIEIs to distributions and plot on top of histograms
            if ~isempty(logIEIsThisCell) && size(logIEIsThisCell, 1) > 1
                [~, ~, ~, ~, params] = ...
                    fit_logIEI(logIEIsThisCell, idCellFull, 'XUnit', xUnitHist, ...
                            'TruncateFlag', truncateDataFlag, ...
                            'PlotFlag', all([plotFitsFlag, plotPerCellFlag]), ...
                            'OutFolder', outFolder);

                % Store parameter fits in data structure
                fitsGrouped.(idGroup).(idSlice).(idCell).logIEIparams = params;

                if isfield(params, 'mu2Model0')
                    ctCellLogModel0 = ctCellLogModel0 + 1;
                    sliceNumLogModel0(ctCellLogModel0, 1) = iSlice;
                    mu1LogModel0AllCells(ctCellLogModel0, 1) = params.mu1Model0;
                    mu2LogModel0AllCells(ctCellLogModel0, 1) = params.mu2Model0;
                    spacingLogModel0AllCells(ctCellLogModel0, 1) = params.spacingModel0;
                    threshold1LogModel0AllCells(ctCellLogModel0, 1) = params.thresholdModel0;
                    threshold2LogModel0AllCells(ctCellLogModel0, 1) = params.thresholdModel0;
                    void1LogModel0AllCells(ctCellLogModel0, 1) = params.voidModel0;
                    void2LogModel0AllCells(ctCellLogModel0, 1) = params.voidModel0;
                    if plotOverlapFlag || plotHeatMapFlag
                        pdfLogModel0AllCells{iSlice} = [pdfLogModel0AllCells{iSlice}; ...
                                                        params.pdfModel0(xVec)];
                    end
                    fprintf('Done fitting %s with a kernel distribution! \n\n', idCellFull);
                else
                    fprintf('%s was not fitted with a kernel distribution! \n\n', idCellFull);
                end
                if isfield(params, 'mu1Model2')
                    ctCellLogModel2 = ctCellLogModel2 + 1;
                    sliceNumLogModel2(ctCellLogModel2, 1) = iSlice;
                    mu1LogModel2AllCells(ctCellLogModel2, 1) = params.mu1Model2;
                    mu2LogModel2AllCells(ctCellLogModel2, 1) = params.mu2Model2;
                    spacingLogModel2AllCells(ctCellLogModel2, 1) = params.spacingModel2;
                    threshold1LogModel2AllCells(ctCellLogModel2, 1) = params.threshold1Model2;
                    threshold2LogModel2AllCells(ctCellLogModel2, 1) = params.threshold2Model2;
                    void1LogModel2AllCells(ctCellLogModel2, 1) = params.void1Model2;
                    void2LogModel2AllCells(ctCellLogModel2, 1) = params.void2Model2;
                    if plotOverlapFlag || plotHeatMapFlag
                        pdfLogModel2AllCells{iSlice} = [pdfLogModel2AllCells{iSlice}; ...
                                                    params.pdfModel2(xVec)];
                        comp2LogModel2AllCells{iSlice} = [comp2LogModel2AllCells{iSlice}; ...
                                                    params.comp2Model2(xVec)];
                    end
                    fprintf('Done fitting %s with two Gaussians! \n\n', idCellFull);
                else
                    fprintf('%s cannot be fitted with two Gaussians! \n\n', ...
                            idCellFull);
                end

                if isfield(params, 'mu1Model3')
                    ctCellLogModel3 = ctCellLogModel3 + 1;
                    sliceNumModel3(ctCellLogModel3, 1) = iSlice;
                    mu1LogModel3AllCells(ctCellLogModel3, 1) = params.mu1Model3;
                    mu2LogModel3AllCells(ctCellLogModel3, 1) = params.mu2Model3;
                    spacingLogModel3AllCells(ctCellLogModel3, 1) = params.spacingModel3;
                    threshold1LogModel3AllCells(ctCellLogModel3, 1) = params.threshold1Model3;
                    threshold2LogModel3AllCells(ctCellLogModel3, 1) = params.threshold2Model3;
                    void1LogModel3AllCells(ctCellLogModel3, 1) = params.void1Model3;
                    void2LogModel3AllCells(ctCellLogModel3, 1) = params.void2Model3;
                    if plotOverlapFlag || plotHeatMapFlag
                        pdfLogModel3AllCells{iSlice} = [pdfLogModel3AllCells{iSlice}; ...
                                                    params.pdfModel3(xVec)];
                        comp2LogModel3AllCells{iSlice} = [comp2LogModel3AllCells{iSlice}; ...
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
            if ~isempty(pdfLogModel0AllCells)
                h = figure('Visible', 'off');
                hold on;
                plot(xVec, pdfLogModel0AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Kernel distributions for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(figsFolder, idGroup, ...
                        [idGroup, '_', idSlice, '_pdfs_Kernel']), 'png');
                close(h);
            end

            if ~isempty(pdfLogModel2AllCells)
                h = figure('Visible', 'off');
                hold on;
                plot(xVec, pdfLogModel2AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Two-Gaussian fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(figsFolder, idGroup, ...
                        [idGroup, '_', idSlice, '_pdfs_GaussOnly']), 'png');
                close(h);

                h = figure('Visible', 'off');
                hold on;
                plot(xVec, comp2LogModel2AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['2nd Gaussian fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(figsFolder, idGroup, ...
                        [idGroup, '_', idSlice, '_comp2s_GaussOnly']), 'png');
                close(h);
            end

            if ~isempty(pdfLogModel3AllCells)
                h = figure('Visible', 'off');
                hold on;
                plot(xVec, pdfLogModel3AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Gaussian-Exp-Exponential fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(figsFolder, idGroup, ...
                        [idGroup, '_', idSlice, '_pdfs_GaussExp']), 'png');
                close(h);

                h = figure('Visible', 'off');
                hold on;
                plot(xVec, comp2LogModel3AllCells);
                xlabel(xUnitHist);
                ylabel('Probability density');
                title(['Exp-Exponential fits for cells of ', idSlice]);
                legend(idAllCells{iSlice}, 'Interpreter', 'none', 'location', 'northwest');
                saveas(h, fullfile(figsFolder, idGroup, ...
                        [idGroup, '_', idSlice, '_comp2s_GaussExp']), 'png');
                close(h);
            end
        end
    end
    % Combine stats
    % commands = generate_stats_commands('Cell', 'Model0');
    % cellfun(@eval, commands);
    % commands = generate_stats_commands('Cell', 'Model2');
    % cellfun(@eval, commands);
    % commands = generate_stats_commands('Cell', 'Model3');
    % cellfun(@eval, commands);
    % commands = generate_stats_commands('Cell', 'Model4');
    % cellfun(@eval, commands);
    commands = generate_stats_commands('Cell', 'LogModel0');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Cell', 'LogModel2');
    cellfun(@eval, commands);
    commands = generate_stats_commands('Cell', 'LogModel3');
    cellfun(@eval, commands);

    for iSlice = 1:nSlices      % for each slice
        % Extract identifier for slice
        idSlice = allSliceLabels{iSlice};
        fprintf('Slice %s ... \n', idSlice);

        % Extract IEIs for each cell of this slice and remove pooled logIEIs
        thisSlice = rmfield(thisGroup.(idSlice), 'interEventIntervals');
                                            % a structure

        % Get all cells
        allCellLabels = fieldnames(thisSlice);
        nCells = numel(allCellLabels);
        
        if nCells > 0 && plotHeatMapFlag
            commands = generate_heatmap_commands('Cell', 'LogModel0', figsFolder);
            cellfun(@eval, commands);
            commands = generate_heatmap_commands('Cell', 'LogModel2', figsFolder);
            cellfun(@eval, commands);
            commands = generate_heatmap_commands('Cell', 'LogModel3', figsFolder);
            cellfun(@eval, commands);
        end
    end

    % Plot histograms of statistics of all cells grouped by slice
    if plotStatHistFlag
        figname = fullfile(figsFolder, ['mu1_', idGroup, '_Kernel']);
        plot_grouped_histogram(figname, exp(mu1LogModel0AllCells), ...
                    sliceNumLogModel0, allSliceLabels, ...
                    'Intra-burst means', 's', ...
                    ['Intra-burst means of ', idGroup, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['mu1_', idGroup, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(mu1LogModel2AllCells), ...
                    sliceNumLogModel2, allSliceLabels, ...
                    'Intra-burst means', 's', ...
                    ['Intra-burst means of ', idGroup, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['mu1_', idGroup, '_GaussExp']);
        plot_grouped_histogram(figname, exp(mu1LogModel3AllCells), ...
                    sliceNumLogModel3, allSliceLabels, ...
                    'Intra-burst means', 's', ...
                    ['Intra-burst means of ', idGroup, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['mu2_', idGroup, '_Kernel']);
        plot_grouped_histogram(figname, exp(mu2LogModel0AllCells), ...
                    sliceNumLogModel0, allSliceLabels, ...
                    'Inter-burst means', 's', ...
                    ['Inter-burst means of ', idGroup, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['mu2_', idGroup, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(mu2LogModel2AllCells), ...
                    sliceNumLogModel2, allSliceLabels, ...
                    'Inter-burst means', 's', ...
                    ['Inter-burst means of ', idGroup, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['mu2_', idGroup, '_GaussExp']);
        plot_grouped_histogram(figname, mu2LogModel3AllCells, ...
                    sliceNumLogModel3, allSliceLabels, ...
                    'Inter-burst means', 's', ...
                    ['Inter-burst means of ', idGroup, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['spacing_', idGroup, '_Kernel']);
        plot_grouped_histogram(figname, spacingLogModel0AllCells, ...
                    sliceNumLogModel0, allSliceLabels, ...
                    'Spacing parameters', 'log s', ...
                    ['Spacing parameters of ', idGroup, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['spacing_', idGroup, '_GaussOnly']);
        plot_grouped_histogram(figname, spacingLogModel2AllCells, ...
                    sliceNumLogModel2, allSliceLabels, ...
                    'Spacing parameters', 'log s', ...
                    ['Spacing parameters of ', idGroup, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['spacing_', idGroup, '_GaussExp']);
        plot_grouped_histogram(figname, spacingLogModel3AllCells, ...
                    sliceNumLogModel3, allSliceLabels, ...
                    'Spacing parameters', 'log s', ...
                    ['Spacing parameters of ', idGroup, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['threshold_', idGroup, '_Kernel']);
        plot_grouped_histogram(figname, exp(threshold1LogModel0AllCells), ...
                    sliceNumLogModel0, allSliceLabels, ...
                    'Thresholds', 's', ...
                    ['Thresholds of ', idGroup, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['threshold1_', idGroup, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(threshold1LogModel2AllCells), ...
                    sliceNumLogModel2, allSliceLabels, ...
                    'Threshold #1s', 's', ...
                    ['Threshold #1s of ', idGroup, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['threshold1_', idGroup, '_GaussExp']);
        plot_grouped_histogram(figname, exp(threshold1LogModel3AllCells), ...
                    sliceNumLogModel3, allSliceLabels, ...
                    'Threshold #1s', 's', ...
                    ['Threshold #1s of ', idGroup, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['threshold2_', idGroup, '_GaussOnly']);
        plot_grouped_histogram(figname, exp(threshold2LogModel2AllCells), ...
                    sliceNumLogModel2, allSliceLabels, ...
                    'Threshold #2s', 's', ...
                    ['Threshold #2s of ', idGroup, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['threshold2_', idGroup, '_GaussExp']);
        plot_grouped_histogram(figname, exp(threshold2LogModel3AllCells), ...
                    sliceNumLogModel3, allSliceLabels, ...
                    'Threshold #2s', 's', ...
                    ['Threshold #2s of ', idGroup, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count');
        figname = fullfile(figsFolder, ['void_', idGroup, '_Kernel']);
        plot_grouped_histogram(figname, void1LogModel0AllCells, ...
                    sliceNumLogModel0, allSliceLabels, ...
                    'Void parameters', '', ...
                    ['Void parameters of ', idGroup, ...
                    ' using kernel distribution fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(figsFolder, ['void1_', idGroup, '_GaussOnly']);
        plot_grouped_histogram(figname, void1LogModel2AllCells, ...
                    sliceNumLogModel2, allSliceLabels, ...
                    'Void parameter #1s', '', ...
                    ['Void parameter #1s of ', idGroup, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(figsFolder, ['void1_', idGroup, '_GaussExp']);
        plot_grouped_histogram(figname, void1LogModel3AllCells, ...
                    sliceNumLogModel3, allSliceLabels, ...
                    'Void parameter #1s', '', ...
                    ['Void parameter #1s of ', idGroup, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(figsFolder, ['void2_', idGroup, '_GaussOnly']);
        plot_grouped_histogram(figname, void2LogModel2AllCells, ...
                    sliceNumLogModel2, allSliceLabels, ...
                    'Void parameter #2s', '', ...
                    ['Void parameter #2s of ', idGroup, ...
                    ' using Two-Gaussian fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
        figname = fullfile(figsFolder, ['void2_', idGroup, '_GaussExp']);
        plot_grouped_histogram(figname, void2LogModel3AllCells, ...
                    sliceNumLogModel3, allSliceLabels, ...
                    'Void parameter #2s', '', ...
                    ['Void parameter #2s of ', idGroup, ...
                    ' using Gaussian-ExpExponential fits'], ...
                    'YLabel', 'Cell count', 'XLimits', [0, 1]);
    end

    % Plot scatter plots of statistics of all cells grouped by slice
    if plotScatterFlag
        commands = generate_scatter_commands(idGroup, 'LogModel0', figsFolder);
        cellfun(@eval, commands);
        commands = generate_scatter_commands(idGroup, 'LogModel2', figsFolder);
        cellfun(@eval, commands);
        commands = generate_scatter_commands(idGroup, 'LogModel3', figsFolder);
        cellfun(@eval, commands);
    end
end

% Plot bar graphs comparing means
if plotBarFlag
    h = figure(1017);
    clf(h);
    h = bar_w_CI(h, muPerCellLogModel0, muPerCellLogModel0Low, ...
                    muPerCellLogModel0High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using kernel distributions for each cell');
    saveas(h, fullfile(figsFolder, 'Kernel_means_AvgCells'), 'png');
    
    h = figure(1018);
    clf(h);
    h = bar_w_CI(h, muPerCellLogModel2, muPerCellLogModel2Low, ...
                    muPerCellLogModel2High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using two-Gaussian fits for each cell');
    saveas(h, fullfile(figsFolder, 'GaussOnly_means_AvgCells'), 'png');

    h = figure(1019);
    clf(h);
    h = bar_w_CI(h, muPerCellLogModel3, muPerCellLogModel3Low, ...
                    muPerCellLogModel3High, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Intra-burst', 'Inter-burst'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'northwest');
    ylabel('Mean inter-event interval (s)');
    title('Mean IEIs using Gaussian-Exp-Exponential fits for each cell');
    saveas(h, fullfile(figsFolder, 'GaussExp_means_AvgCells'), 'png');

    h = figure(1020);
    clf(h);
    h = bar_w_CI(h, [muPerCellLogModel0(1, :); muPerCellLogModel2(1, :); ...
                        muPerCellLogModel3(1, :)], ...
                    [muPerCellLogModel0Low(1, :); muPerCellLogModel2Low(1, :); ...
                        muPerCellLogModel3Low(1, :)], ...
                    [muPerCellLogModel0High(1, :); muPerCellLogModel2High(1, :); ...
                        muPerCellLogModel3High(1, :)], ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Mean inter-spike interval (s)');
    title('Mean intraburst ISIs averaged across cells');
    saveas(h, fullfile(figsFolder, 'Mu1_AvgCells'), 'png');

    h = figure(1021);
    clf(h);
    h = bar_w_CI(h, spacingPerCell, spacingPerCellLow, spacingPerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Time scale spacing parameter');
    title('Spacing parameters averaged across cells');
    saveas(h, fullfile(figsFolder, 'Spacing_AvgCells'), 'png');

    h = figure(1022);
    clf(h);
    h = bar_w_CI(h, threshold1PerCell, threshold1PerCellLow, threshold1PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #1 (s)');
    title('Burst detection threshold #1 averaged across cells');
    saveas(h, fullfile(figsFolder, 'Threshold1_AvgCells'), 'png');

    h = figure(1023);
    clf(h);
    h = bar_w_CI(h, threshold2PerCell, threshold2PerCellLow, threshold2PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Burst detection threshold #2 (s)');
    title('Burst detection threshold #2 averaged across cells');
    saveas(h, fullfile(figsFolder, 'Threshold2_AvgCells'), 'png');

    h = figure(1024);
    clf(h);
    h = bar_w_CI(h, void1PerCell, void1PerCellLow, void1PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #1');
    title('Void parameter #1s averaged across cells');
    saveas(h, fullfile(figsFolder, 'Void1_AvgCells'), 'png');

    h = figure(1025);
    clf(h);
    h = bar_w_CI(h, void2PerCell, void2PerCellLow, void2PerCellHigh, ...
                    'BarSeparation', 0.135, ...
                    'CILineWidth', 2, ...
                    'XTickLabel', {'Kernel', 'Gauss-Gauss', 'Gauss-ExpExp'});
    legend(allGroups, 'Interpreter', 'none', 'location', 'eastoutside');
    ylabel('Void parameter #2');
    title('Void parameter #2s averaged across cells');
    saveas(h, fullfile(figsFolder, 'Void2_AvgCells'), 'png');
end

% Save fitsGrouped structure as a matfile
save(fullfile(outFolderParent, [outFileBase, '.mat']), ...
        '-struct', 'fitsGrouped');

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
    idMethod = 'GaussExp';
    iModel = 3;
case 'Model4'
    idMethod = 'Bodova';
    iModel = 4;
case 'LogModel0'
    idMethod = 'LogKernel';
    iModel = 5;
case 'LogModel2'
    idMethod = 'LogGaussOnly';
    iModel = 6;
case 'LogModel3'
    idMethod = 'LogGaussExpExp';
    iModel = 7;
otherwise
    error('ID of model not recognized!');
end

commands{1} = sprintf([...
    'meanMu1', idModel, ' = nanmean(mu1', idModel, 'All', perWhat, 's);\n', ...
    'meanMu2', idModel, ' = nanmean(mu2', idModel, 'All', perWhat, 's);\n', ...
    'stdErrMu1', idModel, ' = nanstderr(mu1', idModel, 'All', perWhat, 's);\n', ...
    'stdErrMu2', idModel, ' = nanstderr(mu2', idModel, 'All', perWhat, 's);\n', ...
    'muPer', perWhat, idModel, '(1, iGroup) = exp(meanMu1', idModel, ');\n', ...
    'muPer', perWhat, idModel, 'Low(1, iGroup) = exp(meanMu1', idModel, ' - 1.96*stdErrMu1', idModel, ');\n', ...
    'muPer', perWhat, idModel, 'High(1, iGroup) = exp(meanMu1', idModel, ' + 1.96*stdErrMu1', idModel, ');\n', ...
    ]);
switch idModel
case {'LogModel0', 'LogModel2'}
    commands{2} = sprintf([...
        'muPer', perWhat, idModel, '(2, iGroup) =  ', ...
            'exp(meanMu2', idModel, ');\n', ...
        'muPer', perWhat, idModel, 'Low(2, iGroup) =  ', ...
            'exp(meanMu2', idModel, ' - 1.96*stdErrMu2', idModel, ');\n', ...
        'muPer', perWhat, idModel, 'High(2, iGroup) =  ', ...
            'exp(meanMu2', idModel, ' + 1.96*stdErrMu2', idModel, ');\n', ...
        ]);
case 'LogModel3'
    commands{2} = sprintf([...
        'muPer', perWhat, idModel, '(2, iGroup) =  ', ...
            'meanMu2', idModel, ';\n', ...
        'muPer', perWhat, idModel, 'Low(2, iGroup) =  ', ...
            'meanMu2', idModel, ' - 1.96*stdErrMu2', idModel, ';\n', ...
        'muPer', perWhat, idModel, 'High(2, iGroup) =  ', ...
            'meanMu2', idModel, ' + 1.96*stdErrMu2', idModel, ';\n', ...
        ]);
otherwise
    error('ID of model not recognized!');
end
commands{3} = sprintf([...
    'meanSpacing', idModel, ' = nanmean(spacing', idModel, 'All', perWhat, 's);\n', ...
    'stdErrSpacing', idModel, ' = nanstderr(spacing', idModel, 'All', perWhat, 's);\n', ...
    'spacingPer', perWhat, '(', num2str(iModel), ', iGroup) = ', ...
		'meanSpacing', idModel, ';\n', ...
    'spacingPer', perWhat, 'Low(', num2str(iModel), ', iGroup) = ', ...
		'meanSpacing', idModel, ' - 1.96*stdErrSpacing', idModel, ';\n', ...
    'spacingPer', perWhat, 'High(', num2str(iModel), ', iGroup) = ', ...
		'meanSpacing', idModel, ' + 1.96*stdErrSpacing', idModel, ';\n', ...
    'meanThreshold1', idModel, ' = nanmean(threshold1', idModel, 'All', perWhat, 's);\n', ...
    'stdErrThreshold1', idModel, ' = nanstderr(threshold1', idModel, 'All', perWhat, 's);\n', ...
    'threshold1Per', perWhat, '(', num2str(iModel), ', iGroup) = ', ...
		'exp(meanThreshold1', idModel, ');\n', ...
    'threshold1Per', perWhat, 'Low(', num2str(iModel), ', iGroup) = ', ...
		'exp(meanThreshold1', idModel, ' - 1.96*stdErrThreshold1', idModel, ');\n', ...
    'threshold1Per', perWhat, 'High(', num2str(iModel), ', iGroup) = ', ...
		'exp(meanThreshold1', idModel, ' + 1.96*stdErrThreshold1', idModel, ');\n', ...
    'meanThreshold2', idModel, ' = nanmean(threshold2', idModel, 'All', perWhat, 's);\n', ...
    'stdErrThreshold2', idModel, ' = nanstderr(threshold2', idModel, 'All', perWhat, 's);\n', ...
    'threshold2Per', perWhat, '(', num2str(iModel), ', iGroup) = ', ...
		'exp(meanThreshold2', idModel, ');\n', ...
    'threshold2Per', perWhat, 'Low(', num2str(iModel), ', iGroup) = ', ...
		'exp(meanThreshold2', idModel, ' - 1.96*stdErrThreshold2', idModel, ');\n', ...
    'threshold2Per', perWhat, 'High(', num2str(iModel), ', iGroup) = ', ...
		'exp(meanThreshold2', idModel, ' + 1.96*stdErrThreshold2', idModel, ');\n', ...
    'meanVoid1', idModel, ' = nanmean(void1', idModel, 'All', perWhat, 's);\n', ...
    'stdErrVoid1', idModel, ' = nanstderr(void1', idModel, 'All', perWhat, 's);\n', ...
    'void1Per', perWhat, '(', num2str(iModel), ', iGroup) = ', ...
		'meanVoid1', idModel, ';\n', ...
    'void1Per', perWhat, 'Low(', num2str(iModel), ', iGroup) = ', ...
		'meanVoid1', idModel, ' - 1.96*stdErrVoid1', idModel, ';\n', ...
    'void1Per', perWhat, 'High(', num2str(iModel), ', iGroup) = ', ...
		'meanVoid1', idModel, ' + 1.96*stdErrVoid1', idModel, ';\n', ...
    'meanVoid2', idModel, ' = nanmean(void2', idModel, 'All', perWhat, 's);\n', ...
    'stdErrVoid2', idModel, ' = nanstderr(void2', idModel, 'All', perWhat, 's);\n', ...
    'void2Per', perWhat, '(', num2str(iModel), ', iGroup) = ', ...
		'meanVoid2', idModel, ';\n', ...
    'void2Per', perWhat, 'Low(', num2str(iModel), ', iGroup) = ', ...
		'meanVoid2', idModel, ' - 1.96*stdErrVoid2', idModel, ';\n', ...
    'void2Per', perWhat, 'High(', num2str(iModel), ', iGroup) = ', ...
		'meanVoid2', idModel, ' + 1.96*stdErrVoid2', idModel, ';\n', ...
    ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function commands = generate_heatmap_commands(perWhat, idModel, figsFolder);

switch idModel
case 'Model0'
    idMethod = 'Kernel';
    descriptionMethod = 'Kernel distributions';
case 'Model2'
    idMethod = 'GaussOnly';
    descriptionMethod = 'Two-Gaussian distributions';
case 'Model3'
    idMethod = 'GaussExp';
    descriptionMethod = 'Gaussian-Exponential distributions';
case 'Model4'
    idMethod = 'Bodova';
    descriptionMethod = 'Bodova et al';
case 'LogModel0'
    idMethod = 'LogKernel';
    descriptionMethod = 'Kernel distributions';
case 'LogModel2'
    idMethod = 'LogGaussOnly';
    descriptionMethod = 'Two-Gaussian distributions';
case 'LogModel3'
    idMethod = 'LogGaussExpExp';
    descriptionMethod = 'Gaussian-ExpExponential distributions';
otherwise
    error('ID of model not recognized!');
end

switch perWhat
case 'Slice'
    str1 = '';
    str2 = '';
    str3 = '';
    str4 = 'idGroup';
case 'Cell'
    str1 = '{iSlice}';
    str2 = ['(sliceNum', idModel, ' == iSlice)'];
    str3 = 'idGroup, ''/'', ';
    str4 = 'idGroup, ''_'', idSlice';
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
        'saveas(h, fullfile(''', figsFolder, ''', [', str3, '''heatmap_'', ', str4, ...
        ', ''_pdfs_', idMethod, ''']), ''png'');\n', ...
        'close(h);\n', ...
    'end', ...
    ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function commands = generate_scatter_commands(idGroup, idModel, figsFolder)

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

figname = fullfile(figsFolder, ['scatter_mu1_mu2_', idGroup, '_', idMethod, '_unscaled']);
commands{1} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(mu1', idModel, 'AllCells), ', str1, idModel, str2, ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Intra-burst mean'', ''s'', ', ...
            '''Inter-burst mean'', ''s'', ', ...
            '''Means for ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false);'];

figname = fullfile(figsFolder, ['scatter_mu1_mu2_', idGroup, '_', idMethod, '_log']);
commands{2} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(mu1', idModel, 'AllCells), ', str1, idModel, str2, ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Intra-burst mean'', ''s'', ', ...
            '''Inter-burst mean'', ''s'', ', ...
            '''Means for ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''log'', ', ...
            '''XLimits'', [0.5, 8], ''YLimits'', [1, 1000]);'];

figname = fullfile(figsFolder, ['scatter_mu1_mu2_', idGroup, '_', idMethod, '_linear']);
commands{3} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(mu1', idModel, 'AllCells), ', str1, idModel, str2, ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Intra-burst mean'', ''s'', ', ...
            '''Inter-burst mean'', ''s'', ', ...
            '''Means for ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false, ', ...
            '''XLimits'', [0.5, 7.5], ''YLimits'', [0, 220]);'];

figname = fullfile(figsFolder, ['scatter_threshold1_threshold2_', idGroup, '_', idMethod, '_unscaled']);
commands{4} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold1', idModel, 'AllCells), exp(threshold2', idModel, 'AllCells), ', ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Threshold (intersection)'', ''s'', ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Thresholds for ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false);'];

figname = fullfile(figsFolder, ['scatter_threshold1_threshold2_', idGroup, '_', idMethod, '_log']);
commands{5} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold1', idModel, 'AllCells), exp(threshold2', idModel, 'AllCells), ', ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Threshold (intersection)'', ''s'', ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Thresholds for ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''log'', ', ...
            '''XLimits'', [2, 70], ''YLimits'', [2, 200]);'];

figname = fullfile(figsFolder, ['scatter_threshold1_threshold2_', idGroup, '_', idMethod, '_linear']);
commands{6} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold1', idModel, 'AllCells), exp(threshold2', idModel, 'AllCells), ', ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Threshold (intersection)'', ''s'', ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Thresholds for ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false, ', ...
            '''XLimits'', [0, 30], ''YLimits'', [0, 40]);'];

figname = fullfile(figsFolder, ['scatter_void1_void2_', idGroup, '_', idMethod]);
commands{7} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'void1', idModel, 'AllCells, void2', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Void Parameter (intersection)'', '''', ', ...
            '''Void Parameter (minimum)'', '''', ', ...
            '''Void Parameters for ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', false, ', ...
            '''XLimits'', [0, 1], ''YLimits'', [0, 1]);'];

figname = fullfile(figsFolder, ['scatter_threshold2_void2_', idGroup, '_', idMethod]);
commands{8} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold2', idModel, 'AllCells), void2', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Void Parameter (minimum)'', '''', ', ...
            '''Void Parameter vs. Threshold for cells in ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''linear'', ', ...
            '''XLimits'', [2, 200], ''YLimits'', [0, 1]);'];

figname = fullfile(figsFolder, ['scatter_threshold2_spacing_', idGroup, '_', idMethod]);
commands{9} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'exp(threshold2', idModel, 'AllCells), spacing', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Threshold (minimum)'', ''s'', ', ...
            '''Spacing Parameter'', '''', ', ...
            '''Spacing Parameter vs. Threshold for cells in ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XScale'', ''log'', ''YScale'', ''linear'', ', ...
            '''XLimits'', [2, 200], ''YLimits'', [0, 5]);'];

figname = fullfile(figsFolder, ['scatter_spacing_void2_', idGroup, '_', idMethod]);
commands{10} = ['plot_grouped_scatter(''', figname, ''', ', ...
            'spacing', idModel, 'AllCells, void2', idModel, 'AllCells, ', ...
            'sliceNum', idModel, ', allSliceLabels, ', ...
            '''Spacing Parameter'', '''', ', ...
            '''Void Parameter (minimum)'', '''', ', ...
            '''Burstiness of cells in ', idGroup, ...
            ' using ', descriptionMethod, ''', ', ...
            '''MarkerLineWidth'', 1.5, ''PlotEllipse'', true, ', ...
            '''XLimits'', [0, 5], ''YLimits'', [0, 1]);'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

IEIs = IEIs(IEIs > 0);

fprintf('logIEIs extracted from %s ... \n', idGroup);

% Mark's Note: I was debugging and noticed that it crashed on line 627 when there was only 1 logIEI (e.g.: x10HC067_3nMAng_data459_Cell12).
% I amended the code so that it skips logIEIs that are EITHER empty OR only
% contain one logIEI (see line 627)
if strcmp('x10HC067_3nMAng_data459_Cell12', idCellFull) == 1
    m = 1
end

IEIs = IEIs(IEIs > 0);

fprintf('logIEIs extracted from %s ... \n', idSliceFull);

% Extract IEIs for each cell of this slice and remove pooled logIEIs
thisSlice = thisGroup.(idSlice).perCell;             % a structure
nCells = numel(thisSlice);
IEIs = thisSlice{iCell};
IEIs = IEIs(IEIs > 0);

fprintf('logIEIs extracted from %s ... \n', idCellFull);

% Extract IEIs for each cell of this slice and remove pooled logIEIs
thisSlice = thisGroup.(idSlice).perCell;             % a cell array
nCells = numel(thisSlice);

%}
