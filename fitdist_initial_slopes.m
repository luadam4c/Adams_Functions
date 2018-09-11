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
% Requires:
%       /home/Matlab/Adams_Functions/fit_kernel.m
%       /home/Matlab/Adams_Functions/fit_gaussians_and_refine_threshold.m
%       /home/Matlab/Adams_Functions/remove_outliers.m
%
% Used by:
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m

% File History:
% 2018-09-11 Moved from /media/adamX/m3ha/data_dclamp/initial_slopes.m

%% Hard-coded parameters (must be consistent with dataDclampExtractor.m)
validMethods = {'kernelVoid', 'threeStdMainComponent'};

%% TODO: Fix
outlierMethod = 'fiveStds';         % method for determining outliers
bandwidth2StdRatio = 1/5;           % ratio of kernel bandwidth to 
                                    %   standard deviation
maxNumComponents = 4;               % maximum number of components
                                    %   if multiple Gaussians are used to fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ThresMethod', thresMethodDefault, ...
    @(x) any(validatestring(x, validMethods)));

% Read from the Input Parser
parse(iP, varargin{:});
thresMethod = validatestring(iP.Results.ThresMethod, validMethods);

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
               'OldThr', 0, 'ThrMin', 0, 'OutFolder', histFolder, ...
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
        peak2s(iNSample) = muBest(indRanked(2));

        % Get mu of the third largest component
        peak3s(iNSample) = muBest(indRanked(3));

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