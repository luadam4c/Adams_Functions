function [rmsGauss, indGauss, message, vectorGauss, ...
                    nKeptWindows, keptWindows] = ...
                compute_rms_Gaussian (dataVector, varargin)
%% Compute the root-mean-square level of the "Gaussian part" of a data vector
% Usage: [rmsGauss, indGauss, message, vectorGauss, ...
%                   nKeptWindows, keptWindows] = ...
%               compute_rms_Gaussian (dataVector, varargin)
% Explanation:
%       1. Partition the data vector into windows
%       2. Shift all windows so that the mean is zero
%       3. Test whether each window is Gaussian
%       4. Put the data points from the Gaussian windows (shifted) together and 
%           take the overall standard deviation
%
% Outputs:
%       rmsGauss    - root-mean-square level of the "Gaussian part"
%                   specified as a numeric scalar
%       indGauss    - indices of the data vector that 
%                       were in the "Gaussian part"
%                   specified as a numeric vector
%       message     - a message with the number of windows used
%                   specified as a character array or a cellstr
%       vectorGauss - the "Gaussian part" of a data vector
%                   specified as a numeric vector
%       nKeptWindows    - the number of windows used
%                       specified as a numeric scalar
%       keptWindows     - the windows used as a matrix
%                       specified as a numeric array
%
% Arguments:    
%       dataVector  - data vector to analyze, usually a time series
%                   must be a numeric vector
%       varargin    - 'WindowSize': size in samples for a window
%                   must be a positive integer scalar
%                   default == minimum of 5 samples or 
%                               however many to make 1000 windows
%                   - 'ZSkewnessThres': skewness cutoff for deciding 
%                                       if a window is Gaussian
%                   must be a nonnegative scalar
%                   default == 2
%                   - 'ZExcessKurtosisThres': excess kurtosis cutoff for 
%                                           deciding if a window is Gaussian
%                   must be a nonnegative scalar
%                   default == 2
%                   - 'NoBiasSkewness': whether to turn OFF bias correction
%                                           for skewness
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false (use bias correction)
%                   - 'NoBiasKurtosis': whether to turn OFF bias correction
%                                           for kurtosis
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false (use bias correction)
%                   - 'NoBiasStd': whether to turn OFF bias correction
%                                      for standard deviation
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false (use bias correction)
%
% Used by:
%       /home/Matlab/Kojis_Functions/find_directional_events.m
%       cd/minEASE_detect_gapfree_events.m
%       cd/minEASE_examine_gapfree_events.m
%       cd/minEASE_gui_examine_events.m
%       /home/EEG_gui/EEG_gui.m

% File History:
% ---------- Created by Koji Takahashi
% 2017-05-24 AL - Annotated
% 2017-05-24 AL - Added noBiasKurtosis, noBiasSkewness, noBiasStd
% 2017-05-24 AL - Bias correction is now applied to kurtosis() and skewness()
% 2017-05-24 AL - Renamed rms_calculator() -> rms_Gaussian()
% 2017-05-24 AL - Added nWindowsDefault and set to 1000
% 2017-05-24 AL - Added minWindowSize and set to 5
% 2017-05-24 AL - Added Input Parser Scheme
% 2018-02-16 AL - Now returns a warning message if no window is Gaussian
% 2018-02-16 AL - Now returns indGauss
% 2018-02-16 AL - Added Explanation
% 2018-02-16 AL - Changed skewnessCutoff and excessKurtosisCutoff to 
%                   zSkewnessThres and zExcessKurtosisThres
% 2018-02-16 AL - Change to use z-score thresholds. 
%                   Add skewnessStderr and excessKurtosisStderr
% 2018-05-21 AL - Changed mean() and std() to nanmean() and nanstd()

%% Defaults if windowSize not provided
minWindowSize   = 5;                % minimum window size (samples) allowed
nWindowsDefault = 1000;             % default number of windows

%% Default skewness and excess kurtosis z-score thresholds
zSkewnessThresDefault       = 2;    % default skewness z-score threshold
                                    %   corresponds to a two tailed test 
                                    %   of ~= 0 at alpha ~ 0.05
                                    % see https://brownmath.com/stat/shape.htm
zExcessKurtosisThresDefault = 2;    % default excess kurtosis z-score threshold
                                    %   corresponds to a two tailed test 
                                    %   of ~= 0 at alpha ~ 0.05
                                    % see https://brownmath.com/stat/shape.htm

%% Default weighting schemes for computing 
%       skewness, kurtosis, standard deviation
%   0 - Bias correction (sample statistics)         (N-1 for the case of stdev)
%   1 - No bias correction (population statistics)  (N for the case of stdev)
noBiasSkewnessDefault   = 0;        % whether to apply population skewness
noBiasKurtosisDefault   = 0;        % whether to apply population kurtosis
noBiasStdDefault        = 0;        % whether to apply population stdev

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help compute_rms_Gaussian'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'compute_rms_Gaussian';

% Add required inputs to an Input Parser
addRequired(iP, 'dataVector', ...               % data vector to analyze
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'WindowSize', [], ...          % window size (samples)
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'ZSkewnessThres', zSkewnessThresDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'ZExcessKurtosisThres', zExcessKurtosisThresDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'NoBiasSkewness', noBiasSkewnessDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'NoBiasKurtosis', noBiasKurtosisDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'NoBiasStd', noBiasStdDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, dataVector, varargin{:});
windowSize           = iP.Results.WindowSize;
zSkewnessThres       = iP.Results.ZSkewnessThres;
zExcessKurtosisThres = iP.Results.ZExcessKurtosisThres;
noBiasSkewness       = iP.Results.NoBiasSkewness;
noBiasKurtosis       = iP.Results.NoBiasKurtosis;
noBiasStd            = iP.Results.NoBiasStd;

%% Extract information from dataVector
ndp = length(dataVector);           % number of data points

%% Generate windows
% Compute default windowSize if not provided
if isempty(windowSize)
    % First try to use the default number of windows
    windowSize = floor(ndp / nWindowsDefault);

    % Then change to minimum window size if necessary
    if windowSize < minWindowSize
        windowSize = minWindowSize;
    end
end

% Determine number of windows from windowSize
nWindows = floor(ndp / windowSize);

% Determine the skewnessStderr and excessKurtosisStderr from windowSize
%   See https://brownmath.com/stat/shape.htm OR
%   Cramer, Duncan. 1997. Basic Statistics for Social Research. Routledge. p. 89
n = windowSize;
skewnessStderr = sqrt( ( 6 * n * (n-1) ) / ( (n-2) * (n+1) * (n+3) ) );
excessKurtosisStderr = 2 * skewnessStderr * ...
                        sqrt( ( (n+1) * (n-1) ) / ( (n-3) * (n+5) ) );

% Restrict to the part of dataVector that will be analyzed
ndpToAnalyze = windowSize*nWindows;         % number of data points to analyze
dataToAnalyze = dataVector(1:ndpToAnalyze); % part of data vector to analyze

% Reorganize dataToAnalyze into a matrix with each window on a separate column
dataMatrix = reshape(dataToAnalyze, windowSize, nWindows);

%% Compute statistics
% Compute the means of each window as row vector
meanWindow = nanmean(dataMatrix);       

% Subtract each column (each window) by its mean
dataMatrix = dataMatrix - repmat(meanWindow, windowSize, 1);

% Find the sample skewness of each column (each window)
%   and find the indices of the ones below skewness cutoff
%   Note: Skewness is a measure of the asymmetry of the data around the mean
%               s = E[(X - mu)^3]/sigma^3
%         A Gaussian distribution has a skewness value of 0
% See https://brownmath.com/stat/shape.htm for the formula for bias correction
skewnessArray = skewness(dataMatrix, noBiasSkewness);
keepForSkewness = find(abs(skewnessArray/skewnessStderr) < zSkewnessThres);

% Find the sample excess kurtosis of each column (each window)
%   and find the indices of the ones below excess kurtosis cutoff
%   Note: Kurtosis is a measure of how outlier-prone a distribution is:
%               k = E[(X - mu)^4]/sigma^4
%         A Gaussian distribution has a kurtosis value of 3,
%         so the "excess kurtosis" is kurtosis minus 3
% See https://brownmath.com/stat/shape.htm for the formula for bias correction
excessKurtosisArray = kurtosis(dataMatrix, noBiasKurtosis) - 3;
keepForKurtosis = find(abs(excessKurtosisArray/excessKurtosisStderr) ...
                        < zExcessKurtosisThres);

% Find the windows that satisfy both skewness & excess kurtosis cutoff
indToKeep = intersect(keepForKurtosis, keepForSkewness);
nKeptWindows = length(indToKeep);           % number of windows kept
if nKeptWindows == 0
    % Return a message saying root-mean-square level returned as NaN
    message = {['There is no window satisfying the skewness ', ...
                'and kurtosis thresholds for normality!'], ...
                'Root-mean-square level returned as NaN!'};
    vectorGauss = [];
    rmsGauss = NaN;
    indGauss = [];
else
    % Return a message showing how many windows were used
    if nKeptWindows == 1
        message = sprintf('%d of %d windows was considered Gaussian', ...
                            nKeptWindows, nWindows);
    else
        message = sprintf('%d of %d windows were considered Gaussian', ...
                            nKeptWindows, nWindows);
    end

    % Find all the windows to use as the "Gaussian part"
    keptWindows = dataMatrix(:, indToKeep);

    % Put the kept windows together into a vector 
    %   (i.e., the "Gaussian part of the original data vector) 
    %   and take the root mean square
    vectorGauss = reshape(keptWindows, 1, windowSize * nKeptWindows);
    rmsGauss = nanstd(vectorGauss, noBiasStd);

    % Find the indices in the original vector that were before each window used
    indBeforeWindow = (indToKeep - 1) * windowSize;

    % Find the indices in the original vector that were used
    indGauss = reshape(repmat(indBeforeWindow, windowSize, 1), ...
                        1, nKeptWindows*windowSize) + ...
                repmat(1:windowSize, 1, nKeptWindows);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:
% AL - dataVector could be a row vector, so this will be inflexible
ndp = size(dataVector, 1);                    % number of data points

% AL - this doesn't seem to be used:
sd = std(dataMatrix);

% AL - loops are slower than matrix operations
for i = 1:nWindows
    dataMatrix(:, i) = dataMatrix(:, i) - meanWindow(i);
end

% AL - the z-score threshold is better because the z-score is dependent
%       on the sample size
keepForSkewness = find(abs(skewnessArray) < skewnessCutoff);
keepForKurtosis = find(abs(excessKurtosisArray) < excessKurtosisCutoff);
%% Default skewness and excess kurtosis cutoffs
skewnessCutoffDefault       = 0.2;  % default skewness cutoff
excessKurtosisCutoffDefault = 0.2;  % default excess kurtosis cutoff

meanWindow = mean(dataMatrix);       
rmsGauss = std(vectorGauss, noBiasStd);


%}
