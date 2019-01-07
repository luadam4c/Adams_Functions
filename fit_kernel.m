function [model, pdfModel, peak1Model, peak2Model, thresholdModel, voidModel, spacingModel] = fit_kernel(X, varargin)
%% Fits a kernel distribution to a data vector and determine the two primary peaks, the threshold and the void and spacing parameters
% Usage: [model, pdfModel, peak1Model, peak2Model, thresholdModel, voidModel, spacingModel] = fit_kernel(X, varargin)
% Explanation:
%       Fit a kernel distribution to a data vector 
%           and determine the two primary peaks,
%           the threshold and the void and spacing parameters
%       The first primary peak is the mode of the distribution and
%           the second primary peak is the peak to the right or left
%           that yields the largest void parameter. This follows from the
%           algorithm of Pasquale et al., 2010
% Example(s):
%       [model, pdfModel, peak1Model, peak2Model, ...
%           thresholdModel, voidModel, spacingModel] = ...
%           fit_kernel(data, 'Bandwidth2StdRatio', 1/5);
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%       output2     - TODO: Description of output2
%                   specified as a TODO
% Arguments:    
%       X           - data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'Bandwidth': kernel distribution bandwidth
%                   must be a positive scalar
%                   default == std(X) * bandwidth2StdRatio
%                   - 'Bandwidth2StdRatio': ratio of kernel bandwidth to 
%                                           standard deviation
%                   must be a positive scalar
%                   default == 1/5
%                   - 'Resolution': resolution for determining the second peak
%                   must be a positive scalar
%                   default == range(X) / 1000
%                   - 'Peak2Direction': Directrion of Peak #2 from the mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'right' - positive x direction
%                       'left'  - negative x direction
%                   default == 'right'
% Requires:
%       /home/Matlab/Adams_Functions/plot_pdf.m
%
% Used by:
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m
%
% File History:
% 2018-06-12 Adapted from code in fit_IEI.m
% 2018-07-22 Improved mode determination
% 2018-07-22 Added the parameter Peak2Direction and allowed the second
%               peak to be on the left
% TODO: Extract code to a function find_mode.m with usage:
%       find_mode(pdfModel, 'NPointsToResolve', nPointsToResolve)
% TODO: Extract code to a function find_second_peak.m with usage:
%       find_second_peak(pdfModel, peak1Model, 'Direction', peak2Direction)
% 

%% Hard-coded parameters
validDirections = {'right', 'left'};

%% Default values for optional arguments
bandwidthDefault = [];
bandwidth2StdRatioDefault = 1/5;    % use a fifth of the standard deviation
                                    %   as the bandwidth by default
resolutionDefault = [];             % default resolution
peak2DirectionDefault = 'right';    % default directrion of Peak2 from the mode

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
addRequired(iP, 'X', ...                    % data to distribute among bins
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Bandwidth', bandwidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Bandwidth2StdRatio', bandwidth2StdRatioDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Resolution', resolutionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Peak2Direction', peak2DirectionDefault, ...
    @(x) any(validatestring(x, validDirections)));

% Read from the Input Parser
parse(iP, X, varargin{:});
kernelBandwidth = iP.Results.Bandwidth;
bandwidth2StdRatio = iP.Results.Bandwidth2StdRatio;
resolution = iP.Results.Resolution;
peak2Direction = validatestring(iP.Results.Peak2Direction, validDirections);

% Set dependent argument defaults
if isempty(kernelBandwidth)
    kernelBandwidth = std(X) * bandwidth2StdRatio;
end
if isempty(resolution)
    resolution = range(X) / 1000;
end

%% Fit to a kernel distribution
% Fit data to a kernel distribution using Gaussian kernels
model = fitdist(X, 'Kernel', 'Bandwidth', kernelBandwidth);

% Get the pdf for the fit
pdfModel = @(x) pdf(model, x);

% Get the number of points to resolve
nPointsToResolve = ceil(range(X)/resolution);

% Get scaled pdf values for the fit
[xValues, pdfValues] = plot_pdf(X, 'PDF', pdfModel, 'PlotFlag', false, ...
                                'NPointsToPlot', nPointsToResolve);

%% Find the mode of the kernel distribution
% First find the approximate mode from the substituted values
[~, idxModeApprox] = max(pdfValues);

% Find all peaks from the substituted values
[~, idxPeaks] = findpeaks(pdfValues);

% Compute the number of peaks
nPeaks = length(idxPeaks);

% Find the peak number of the approximate mode
peaknoModeApprox = find(idxPeaks == idxModeApprox, 1);

% Find the index of the previous peak (if exists) or the first index
if peaknoModeApprox > 1
    idxPreviousPeak = idxPeaks(peaknoModeApprox - 1);
else
    idxPreviousPeak = 1;
end

% Find the index of the next peak (if exists) or the last index
if peaknoModeApprox < nPeaks
    idxNextPeak = idxPeaks(peaknoModeApprox + 1);
else
    idxNextPeak = length(xValues);
end

% Find the x value of the minimum right before the mode
[~, idxMinBeforeModeApprox] =  min(pdfValues(idxPreviousPeak:idxModeApprox));
xMinBeforeModeApprox = xValues((idxPreviousPeak - 1) + idxMinBeforeModeApprox);

% Find the x value of the minimum right after the mode
[~, idxMinAfterModeApprox] =  min(pdfValues(idxModeApprox:idxNextPeak));
xMinAfterModeApprox = xValues((idxModeApprox- 1) + idxMinAfterModeApprox);

% Find the exact mode using fminbnd
peak1Model = fminbnd(@(x) -pdfModel(x), ...
                    xMinBeforeModeApprox, xMinAfterModeApprox);

%% Find the second peak and the threshold
% Find the indices in which the second peak will be found
if strcmp(peak2Direction, 'right')
    indToFind = xValues >= peak1Model;
elseif strcmp(peak2Direction, 'left')
    indToFind = xValues <= peak1Model;
else
    error(['The direction', peak2Direction, ' is not recognized!']);
end

% Extract the corresponding x values and pdf values
xToFind = xValues(indToFind);
pdfToFind = pdfValues(indToFind);

% Find all peaks in this range
[pdfPksToFind, xPksToFind] = findpeaks(pdfToFind, xToFind);

% Select the peak that yields the largest void parameter when paired 
%   with the mode, following the algorithm of Pasquale et al., 2010
nPeaks = length(xPksToFind);
if nPeaks == 0
    % If there are no more peaks, there cannot be a peak #2 or a threshold
    peak2Model = NaN;
    thresholdModel = NaN;

    % In this case, the void and spacing parameters are set to a minimum (0)
    voidModel = 0;
    spacingModel = 0;
else
    % Otherwise, compute thresholds and void parameters for all possible 
    %   peak #2s
    peak2s = zeros(nPeaks, 1);         % stores all possible peak #2s
    thresholds = zeros(nPeaks, 1);     % stores all possible thresholds
    voids = zeros(nPeaks, 1);          % stores all void parameters
    for iToFind = 1:nPeaks
        % Get the x value of the current peak
        peak2s(iToFind) = xPksToFind(iToFind);

        % The threshold is the minimum between the mode and the current peak
        if peak1Model <= peak2s(iToFind)
            thresholds(iToFind) = ...
                fminbnd(pdfModel, peak1Model, peak2s(iToFind));
        else
            thresholds(iToFind) = ...
                fminbnd(pdfModel, peak2s(iToFind), peak1Model);
        end

        % Compute the void parameter corresponding to this pair of peaks
        voids(iToFind) = ...
            1 - pdfModel(thresholds(iToFind)) / ...
                sqrt(pdfModel(peak1Model) * pdfPksToFind(iToFind));
    end
    % Select the peak that gives the largest void parameter
    [voidModel, iSelected] = max(voids);   

    % Record the x value and threshold for this peak
    peak2Model = peak2s(iSelected);
    thresholdModel = thresholds(iSelected);    

    % Compute the spacing parameter for this peak
    spacingModel = abs(peak2Model - peak1Model);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the mean, minimum and maximum of the data
meanData = mean(X);
minData = min(X);
maxData = max(X);

% Compute the distance between the approximate mode and the mean
distMode2Mean = abs(modeModelApprox - meanData);

% Look for a valid region within distMode2Mean of the mode
leftBnd = max(minData, modeModelApprox - distMode2Mean);
rightBnd = min(maxData, modeModelApprox + distMode2Mean);

% Find the exact mode using fminbnd
peak1Model = fminbnd(@(x) -pdfModel(x), leftBnd, rightBnd);

[pdfValuesPeaks, xValuesPeaks] = findpeaks(pdfValues, xValues);
nPeaks = length(xValuesPeaks);
[~, peaknoModeApprox] = max(pdfValuesPeaks);
modeModelApprox = xValuesPeaks(peaknoModeApprox);
% Find the x value of the previous and next peaks, 
%   or the first and last x values, respectively
if peaknoModeApprox > 1
    xPreviousPeak = xValuesPeaks(peaknoModeApprox - 1);
else
    xNextPeak = xValues(1);
end
if peaknoModeApprox < nPeaks
    xNextPeak = xValuesPeaks(peaknoModeApprox + 1);
else
    xNextPeak = xValues(end);
end

[pdfValueModeApprox, idxModeApprox] = max(pdfValues);
[pdfValuesPeaks, idxPeaks] = findpeaks(pdfValues);
peaknoModeApprox = find(pdfValuesPeaks == pdfValueModeApprox, 1);
modeModelApprox = xValues(idxModeApprox);

% Find all peaks to the right of the mode of the kernel distribution
toTheRight = xValues >= peak1Model;         % whether the index of x is 
                                            % to the right of the mode
xToTheRight = xValues(toTheRight);          % parts of x that is 
                                            % to the right of the mode
pdfToTheRight = pdfValues(toTheRight);      % parts of the pdf that is 
                                            % to the right of the mode
[pksToTheRight, xPksToTheRight] = findpeaks(pdfToTheRight, xToTheRight);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%