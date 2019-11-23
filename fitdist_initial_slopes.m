function [bestModels, pdfModels, peak1s, peak2s, peak3s, ...
            threshold1s, threshold2s] = ...
            fitdist_initial_slopes(allNSamples, avgSlopes, varargin)
%% Fits initial slope distributions
% Usage: [bestModels, pdfModels, peak1s, peak2s, peak3s, ...
%           threshold1s, threshold2s] = ...
%           fitdist_initial_slopes(allNSamples, avgSlopes, varargin)
% Arguments:
%       TODO
%       varargin    - 'ThresMethod': threshold method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'kernelVoid'    - use a kernel distribution
%                       'threeStdMainComponent' 
%                                       - use the center Gaussian distribution
%                   default == 'threeStdMainComponent'
%                   - 'OutlierMethod': method for determining outliers
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'boxplot'   - same as the Matlab function boxplot()
%                       'isoutlier' - Use the built-in isoutlier function
%                       'fiveStds'  - Take out data points 
%                                       more than 5 standard deviations away
%                       'threeStds' - Take out data points 
%                                       more than 3 standard deviations away
%                       'twoStds'   - Take out data points 
%                                       more than 2 standard deviations away
%                   default == 'isoutlier'
%                   - 'Bandwidth2StdRatio': ratio of kernel bandwidth to 
%                                           standard deviation
%                   must be a positive scalar
%                   default == 1/5
%                   - 'MaxNumComponents': maximum number of Gaussian components
%                   must be a positive integer scalar
%                   default == 3
%                   - 'OutFolder': directory to save figures
%                   must be a string scalar or a character vector
%                   default == pwd
%
%
% Requires:
%       cd/fit_kernel.m
%       cd/fit_gaussians_and_refine_threshold.m
%       cd/remove_outliers.m
%
% Used by:
%       cd/m3ha_initial_slopes.m

% File History:
% 2018-09-11 Moved from cd/m3ha_initial_slopes.m

%% Hard-coded parameters (must be consistent with dataDclampExtractor.m)
validMethods = {'kernelVoid', 'threeStdMainComponent'};
validOutlierMethods = {'boxplot', 'isoutlier', ...
                        'fiveStds', 'threeStds', 'twoStds'};

%% Default values for optional arguments
thresMethodDefault = 'threeStdMainComponent';
outlierMethodDefault = 'fiveStds';  % use five standard deviations from the mean
bandwidth2StdRatioDefault = 1/5;    % use a fifth of the standard deviation
                                    %   as the bandwidth by default
maxNumComponentsDefault = 3;        % maximum number of components
                                    %   if multiple Gaussians are used to fit
outFolderDefault = '';              % default directory to save figures


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ThresMethod', thresMethodDefault, ...
    @(x) any(validatestring(x, validMethods)));
addParameter(iP, 'OutlierMethod', outlierMethodDefault, ...
    @(x) any(validatestring(x, validOutlierMethods)));
addParameter(iP, 'Bandwidth2StdRatio', bandwidth2StdRatioDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'MaxNumComponents', maxNumComponentsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'positive', 'integer'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, varargin{:});
thresMethod = validatestring(iP.Results.ThresMethod, validMethods);
outlierMethod = validatestring(iP.Results.OutlierMethod, validOutlierMethods);
bandwidth2StdRatio = iP.Results.Bandwidth2StdRatio;
maxNumComponents = iP.Results.MaxNumComponents;
outFolder = iP.Results.OutFolder;

% Set dependent argument defaults
if isempty(outFolder)
    outFolder = pwd;
end

% Count the number of possible nSamples to use
nNSamples = length(allNSamples);

% Generate pdfs of distribution fits
bestModels = cell(nNSamples, 1);
pdfModels = cell(nNSamples, 1);
peak1s = zeros(nNSamples, 1);
peak2s = zeros(nNSamples, 1);
peak3s = zeros(nNSamples, 1);
threshold1s = zeros(nNSamples, 1);
threshold2s = zeros(nNSamples, 1);
parfor iNSample = 1:nNSamples
%for iNSample = 1:nNSamples
    % Get the current nSamples
    nSamples = allNSamples(iNSample);

    % Get the current average slopes
    avgSlopesThis = avgSlopes(:, iNSample);

    % Remove outliers if any
    avgSlopesThisTrunc = ...
        remove_outliers(avgSlopesThis, 'OutlierMethod', outlierMethod);

    % Fit to a distribution
    if strcmp(thresMethod, 'kernelVoid')
        % Fit to a kernel distribution and find the threshold on the left side
        [~, pdfModels{iNSample}, peak1s(iNSample), ...
            peak2s(iNSample), threshold1s(iNSample), ~, ~] = ...
            fit_kernel(avgSlopesThisTrunc, 'Peak2Direction', 'Left', ...
                        'Bandwidth2StdRatio', bandwidth2StdRatio);

        % Fit to a kernel distribution and find the threshold on the right side
        [~, pdfModels{iNSample}, peak1s(iNSample), ...
            peak3s(iNSample), threshold2s(iNSample), ~, ~] = ...
            fit_kernel(avgSlopesThisTrunc, 'Peak2Direction', 'Right', ...
                        'Bandwidth2StdRatio', bandwidth2StdRatio);
    elseif strcmp(thresMethod, 'threeStdMainComponent')
        % Fit to a Gaussian distribution

        % Create a figure name for fit_gaussians_and_refine_threshold.m
        figname = ['initial_slope_nSamples_', ...
                        num2str(nSamples), '_', thresMethod, '.png'];

        % Fit Gaussians and refine thresholds
        [bestModel, ~, muBest, ~, ~, ~, ...
            rightThres, leftThres, ~, ~, indRanked] = ...
            fit_gaussians_and_refine_threshold(avgSlopesThisTrunc, figname, ...
               'Slope (V/s)', 'MaxNumComponents', maxNumComponents, ...
               'OldThr', 0, 'ThrMin', 0, 'OutFolder', outFolder, ...
               'FitMode', 0, 'PeakClass', ones(size(avgSlopesThisTrunc)), ...
               'PeakClassLabels', {'data'}, ...
               'MinThreshold', 0, 'ThresMode', thresMethod);

        % Store the best model
        bestModels{iNSample} = bestModel;

        % Get pdf of most representative component
        pdfModels{iNSample} = @(x) pdf(bestModel, x);

        % Get mu of the largest component
        peak1s(iNSample) = muBest(indRanked(1));

        % Get mu of the second largest component
        if length(indRanked) >= 2
            peak2s(iNSample) = muBest(indRanked(2));
        else
            peak2s(iNSample) = NaN;
        end

        % Get mu of the third largest component
        if length(indRanked) >= 3
            peak3s(iNSample) = muBest(indRanked(3));
        else
            peak3s(iNSample) = NaN;
        end

        % Set 1st threshold to 3 standard deviations to the left of the mean
        threshold1s(iNSample) = leftThres;

        % Set 2nd threshold to 3 standard deviations to the right of the mean
        threshold2s(iNSample) = rightThres;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%