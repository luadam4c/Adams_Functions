function [model, pdfModel, peak1Model, peak2Model, thresholdModel, voidModel, spacingModel] = fit_kernel(X, varargin)
%% Fits a kernel distribution to a data vector and determine the two primary peaks, the threshold and the void and spacing parameters
% Usage: [model, pdfModel, peak1Model, peak2Model, thresholdModel, voidModel, spacingModel] = fit_kernel(X, varargin)
% Explanation:
%       Fit a kernel distribution to a data vector 
%           and determine the two primary peaks,
%           the threshold and the void and spacing parameters
%       The first primary peak is the mode of the distribution and
%           the second primary peak is the peak to the right (TODO: make this a parameter)
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
% Requires:
%       /home/Matlab/Adams_Functions/plot_pdf.m
%
% Used by:
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m
%
% File History:
% 2018-06-12 Adapted from code in fit_IEI.m
% 

%% Default values for optional arguments
bandwidthDefault = [];
bandwidth2StdRatioDefault = 1/5;    % use a fifth of the standard deviation
                                    %   as the bandwidth by default
resolutionDefault = [];             % default resolution

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

% Read from the Input Parser
parse(iP, X, varargin{:});
kernelBandwidth = iP.Results.Bandwidth;
bandwidth2StdRatio = iP.Results.Bandwidth2StdRatio;
resolution = iP.Results.Resolution;

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
pdfModel = @(X) pdf(model, X);

% Get the number of points to resolve
nPointsToResolve = ceil(range(X)/resolution);

% Get scaled pdf values for the fit
[xValues, pdfValues] = plot_pdf(X, 'PDF', pdfModel, 'PlotFlag', false, ...
                                'NPointsToPlot', nPointsToResolve);

%% Find the mode of the kernel distribution
% First find the approximate mode from the substituted values
[~, idxModeApprox] = max(pdfValues);
modeModelApprox = xValues(idxModeApprox);

% Get the mean, minimum and maximum of the data
meanData = mean(X);
minData = min(X);
maxData = max(X);

% Compute the distance between the approximate mode and the mean
distMode2Mean = abs(modeModelApprox - meanData);

% Find the exact mode using fminbnd
leftBnd = max(minData, modeModelApprox - distMode2Mean);
rightBnd = min(maxData, modeModelApprox + distMode2Mean);
peak1Model = fminbnd(@(x) -pdfModel(x), leftBnd, rightBnd);

%% Find the second peak and the threshold
% Find all peaks to the right of the mode of the kernel distribution
%   TODO: Make "to the right" a parameter
toTheRight = xValues >= peak1Model;         % whether the index of x is 
                                            % to the right of the mode
xToTheRight = xValues(toTheRight);          % parts of x that is 
                                            % to the right of the mode
pdfToTheRight = pdfValues(toTheRight);      % parts of the pdf that is 
                                            % to the right of the mode
[pksToTheRight, locsToTheRight] = findpeaks(pdfToTheRight);

% Select the peak that gives the largest void parameter 
%   when paired with the mode
%   This follows the algorithm of Pasquale et al., 2010
nToTheRight = length(locsToTheRight);
if nToTheRight == 0
    peak2Model = NaN;
    thresholdModel = NaN;
    voidModel = 0;
    spacingModel = 0;
else
    thresholds = zeros(nToTheRight, 1);     % stores all possible thresholds
    voids = zeros(nToTheRight, 1);          % stores all void parameters
    for iToTheRight = 1:nToTheRight
        % The threshold is the minimum between the mode and the current peak
        thresholds(iToTheRight) = ...
            fminbnd(pdfModel, peak1Model, ...
                    xToTheRight(locsToTheRight(iToTheRight)));

        % Compute the void parameter corresponding to this pair of peaks
        voids(iToTheRight) = ...
            1 - pdfModel(thresholds(iToTheRight)) / ...
                sqrt(pdfModel(peak1Model) * pksToTheRight(iToTheRight));
    end
    % Select the peak that gives the largest void parameter
    [voidModel, iSelected] = max(voids);   

    % Record the x value and threshold for this peak
    peak2Model = xToTheRight(locsToTheRight(iSelected));
    thresholdModel = thresholds(iSelected);    

    % Compute the spacing parameter for this peak
    spacingModel = peak2Model - peak1Model;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}