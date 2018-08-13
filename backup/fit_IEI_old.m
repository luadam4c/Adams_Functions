function [model1, model2, model3, outparams] = fit_IEI (data, identifier, varargin)
%% Fit IEI data to curves
% Usage: [model1, model2, model3, outparams] = fit_IEI (data, identifier, varargin)
% Explanation:
%       TODO
% Side Effects:
%       Plots TODO
% Arguments:
%       data        - a vector of IEIs
%                   must be a numeric nonnegative vector
%       identifier  - an identifier for plotting data
%                   must be a string scalar or a character vector
%       varargin    - 'XUnit': unit for IEIs
%                   must be a string scalar or a character vector
%                   default == ms
%                   - 'XLimits': x limits of histogram
%                   must be a 2-element number vector
%                   default == [min(data) max(data)] 
%                   - 'TruncateFlag': whether to truncate data
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': directory for saving plots
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'PlotFlag': whether to plot histograms and fits
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Used By:
%       /media/adamX/Paula_IEIs/paula_iei3.m
%
% File History: 
% 2017-10-19 - Moved from paula_isi2.m
% 2017-10-19 - Changed cutoff to parameter xLimits
% 2017-10-19 - Added XUnit and TruncateFlag

%% Parameters
nPoints = 500;                      % number of points for plotting pdfs
defaultLineWidth = 2;               % default line width for plots
linewidthLines = 0.5;               % line width for lines
colorHist = 'k';                    % color of histogram
lambdaInit = 0.7;                   % initial weight of Gaussian part
                                    %   in Gaussian-Exponential Fit

%% Default values for optional arguments
xUnitDefault = '';                  % default unit for IEIs
xLimitsDefault = [];                % default x limits of histogram
truncateFlagDefault = false;        % whether to truncate data by default
outFolderDefault = '';              % default directory for saving plots
plotFlagDefault = true;             % whether to plot by default
skipModel3Default = false;          % whether to skip model 3 by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help fit_IEI'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'fit_IEI';

% Add required inputs to the Input Parser
addRequired(iP, 'data', ...                 % a vector of IEIs
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'vector'}));
addRequired(iP, 'identifier', ...           % an identifier for plotting data
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XUnit', xUnitDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'TruncateFlag', truncateFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SkipModel3', skipModel3Default, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, data, identifier, varargin{:});
xUnit        = iP.Results.XUnit;
xLimits      = iP.Results.XLimits;
truncateFlag = iP.Results.TruncateFlag;
outFolder    = iP.Results.OutFolder;
plotFlag     = iP.Results.PlotFlag;
skipModel3   = iP.Results.SkipModel3;

% Set dependent argument defaults
if isempty(xUnit)
    xUnit = 'ms';
    fprintf('Unit for IEIs not provided. Set to default %s!\n', xUnit);
end
if isempty(xLimits)
    xLimits = [min(data), max(data)];   
end
if isempty(outFolder)
    outFolder = pwd;
end

%% Fit data
fprintf('Fitting data for %s ... \n', identifier);
if truncateFlag
    % Eliminate data outside of xLimits
    data = data(data <= xLimits(2));
    data = data(data >= xLimits(1));
end

% Total data count and mean of data
% Calculate statistics from data
npts = length(data);
meanIEI = mean(data);
medianIEI = median(data);

% New IEI sequence
%   Chen et al., 2009 
%   ("Detection of bursts in neuronal spike trains 
%       by the mean inter-event interval method")
lessThanMean = data(data < meanIEI);
meanSmallIEI = mean(lessThanMean);

% Number of bins is square root of npts data count
%   Cocatre-Zilgien & Delcomyn, 1992 
%   ("Identification of bursts in spike trains")
% nBins = round(sqrt(npts));
nBins = 1.87 * (npts - 1) ^ 0.4;

% Compute histogram edges, area, and x values for fit pdf
edges = linspace(xLimits(1), xLimits(2), nBins);
binWidth = xLimits(2)/nBins;        % bin width
harea = binWidth * npts;            % total histogram area
x = linspace(xLimits(1), xLimits(2), nPoints)';

try
    % Fit data to one Gaussian
    model1 = fitgmdist(data, 1);

    % Find mean and standard deviations for the fits
    muModel1 = model1.mu;
    sigmaModel1 = model1.Sigma;

    % Get pdf for the fit
    pdfModel1 = harea * pdf(model1, x);
catch
    model1 = [];
end

try
    % Fit data to two Gaussians
    model2 = fitgmdist(data, 2);

    % Find mean and standard deviations for the fits
    [~, origInd] = sort(model2.mu, 'ascend');
    mu1Model2 = model2.mu(origInd(1));
    mu2Model2 = model2.mu(origInd(2));
    sigma1Model2 = model2.Sigma(origInd(1));
    sigma2Model2 = model2.Sigma(origInd(2));
    prop1Model2 = model2.ComponentProportion(origInd(1));
    prop2Model2 = model2.ComponentProportion(origInd(2));

    % Get pdf for the fit
    pdfModel2 = harea * pdf(model2, x);

    % Get pdfs for components of the mixture fit
    comp1Model2 = prop1Model2 * harea * normpdf(x, mu1Model2, sigma1Model2);
    comp2Model2 = prop2Model2 * harea * normpdf(x, mu2Model2, sigma2Model2);

catch
    model2 = [];
end

try
    % Fit data to a mixture of Gaussian and Exponential
    mypdf = @(X, lambda, mu1, sigma, mu2) ...
                    lambda * normpdf(X, mu1, sigma) + ...
                    (1 - lambda) * exppdf(X, mu2);
    start = [lambdaInit, mu1Model2, sigma1Model2, mu2Model2];
    [phat, pci] = mle(data, 'pdf', mypdf, 'start', start);
    lambdaModel3 = phat(1);
    lambdaCIModel3 = pci(:, 1);
    mu1Model3 = phat(2);
    mu1CIModel3 = pci(:, 2);
    sigmaModel3 = phat(3);
    sigmaCIModel3 = pci(:, 3);
    mu2Model3 = phat(4);
    mu2CIModel3 = pci(:, 4);
    model3 = @(X) mypdf(X, lambdaModel3, mu1Model3, sigmaModel3, mu2Model3);
    pdfModel3 = harea * model3(x);
    comp1Model3 = lambdaModel3 * harea * normpdf(x, mu1Model3, sigmaModel3);
    comp2Model3 = (1 - lambdaModel3) * harea * exppdf(x, mu2Model3);
catch
    model3 = [];
end

%% Plot histograms and fits
if plotFlag
    fprintf('Plotting data for %s ... \n', identifier);

    % Change default values
    set(groot, 'defaultLineLineWidth', defaultLineWidth);

    if ~isempty(model1)
        % Plot histograms with fits
        h = figure(1001);
        clf(h);
        hold on;
        hist(data, edges, 'FaceColor', colorHist);
        xlim(xLimits);
        yLimits = get(gca, 'YLim');
        line(meanIEI * ones(1, 2), yLimits, 'Color', 'g', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mean = ', num2str(meanIEI, 3)]);
        line(medianIEI * ones(1, 2), yLimits, 'Color', 'm', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Median = ', num2str(medianIEI, 3)]);
        line(meanSmallIEI * ones(1, 2), yLimits, 'Color', 'r', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Chen threshold = ', num2str(meanSmallIEI, 3)]);
        legend('location', 'northeast');
        xlabel(['Inter-event intervals (', xUnit, ')']);
        ylabel('Count');
        title(['Raw Data for ', strrep(identifier, '_', '\_')]);
        saveas(h, fullfile(outFolder, [identifier, '_HistOnly']), 'png');
    end

    if ~isempty(model2)
        % Plot histograms with Gaussian fits
        h = figure(1002);
        clf(h);
        hold on;
        hist(data, edges, 'FaceColor', colorHist);
        xlim(xLimits);
        plot(x, pdfModel1, 'g', 'Displayname', 'single Gaussian');
        plot(x, pdfModel2, 'c', 'Displayname', 'double Gaussian');
        plot(x, comp1Model2, 'm', 'Displayname', 'Component #1 of double Gaussian');
        plot(x, comp2Model2, 'r', 'Displayname', 'Component #2 of double Gaussian');
        yLimits = get(gca, 'YLim');
        line(muModel1 * ones(1, 2), yLimits, 'Color', 'g', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mean of Single Gaussian = ', num2str(muModel1, 3)]);
        line(mu1Model2 * ones(1, 2), yLimits, 'Color', 'm', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mean of Component #1 = ', num2str(mu1Model2, 3)]);
        line(mu2Model2 * ones(1, 2), yLimits, 'Color', 'r', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mean of Component #2 = ', num2str(mu2Model2, 3)]);
        legend('location', 'northeast');
        xlabel(['Inter-event intervals (', xUnit, ')']);
        ylabel('Count');
        title(['Gaussian Fits for ', strrep(identifier, '_', '\_')]);
        saveas(h, fullfile(outFolder, [identifier, '_GaussOnly']), 'png');
    end
    
    if ~isempty(model3)
        % Plot histograms with Gaussian-Exponential fits
        h = figure(1003);
        clf(h);
        hold on;
        hist(data, edges, 'FaceColor', colorHist);
        xlim(xLimits);
        plot(x, pdfModel3, 'c', 'Displayname', 'Gaussian + Exponential');
        plot(x, comp1Model3, 'm', 'Displayname', 'Gaussian part');
        plot(x, comp2Model3, 'r', 'Displayname', 'Exponential part');
        yLimits = get(gca, 'YLim');
        line(mu1Model3 * ones(1, 2), yLimits, 'Color', 'm', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mean of Gaussian part = ', num2str(mu1Model3, 3)]);
        line(mu2Model3 * ones(1, 2), yLimits, 'Color', 'r', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mean of Exponential part = ', num2str(mu2Model3, 3)]);
        legend('location', 'northeast');
        xlabel(['Inter-event intervals (', xUnit, ')']);
        ylabel('Count');
        title(['Gaussian-Exponential Fit for ', strrep(identifier, '_', '\_')]);
        saveas(h, fullfile(outFolder, [identifier, '_GaussExp']), 'png');

        % Remove default values
        set(groot, 'defaultLineLineWidth', 'remove');
    end
end

%% Package parameters into outparams
outparams.meanIEI = meanIEI;
outparams.medianIEI = medianIEI;
outparams.meanSmallIEI = meanSmallIEI;

if ~isempty(model1)
    outparams.pdfModel1 = pdfModel1;
    outparams.muModel1 = muModel1;
    outparams.sigmaModel1 = sigmaModel1;
end

if ~isempty(model2)
    outparams.pdfModel2 = pdfModel2;
    outparams.mu1Model2 = mu1Model2;
    outparams.mu2Model2 = mu2Model2;
    outparams.sigma1Model2 = sigma1Model2;
    outparams.sigma2Model2 = sigma2Model2;
    outparams.prop1Model2 = prop1Model2;
    outparams.prop2Model2 = prop2Model2;
end

if ~isempty(model3)
    outparams.pdfModel3 = pdfModel3;
    outparams.lambdaModel3 = lambdaModel3;
    outparams.lambdaCIModel3 = lambdaCIModel3;
    outparams.mu1Model3 = mu1Model3;
    outparams.mu1CIModel3 = mu1CIModel3;
    outparams.sigmaModel3 = sigmaModel3;
    outparams.sigmaCIModel3 = sigmaCIModel3;
    outparams.mu2Model3 = mu2Model3;
    outparams.mu2CIModel3 = mu2CIModel3;
    outparams.comp1Model3 = comp1Model3;
    outparams.comp2Model3 = comp2Model3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
    if ~skipModel3
    end
    ~skipModel3 && 
%}
