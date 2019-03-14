function [newData, rowsToKeep] = remove_outliers (oldData, varargin)
%% Removes outliers from a data matrix and return a new matrix
% Usage: [newData, rowsToKeep] = remove_outliers (oldData, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       X = magic(5); remove_outliers(X); X(1) = 100; remove_outliers(X)
%
% Outputs:
%       newData     - data matrix with outlying data points removed
%                   specified as a numeric array
%       rowsToKeep  - row indices of original data matrix to keep
%                   specified as a positive integer array
% Arguments:    
%       oldData     - a data matrix with each column being a condition 
%                       and each row being a data point
%                   must be a numeric, logical, datetime or duration array
%       varargin    - 'WL2IQR': the ratio of whisker length to 
%                                   interquartile range
%                   must be a numeric positive scalar
%                   default == 1.5 (same as the Matlab function boxplot())
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
%                   - 'ReplaceWithNans': whether to replace with NaNs
%                                           instead of removal
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotFlag': whether to plot box plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/force_column_vector.m
%
% Used by:
%       cd/compute_stats.m
%       cd/plot_histogram.m
%       cd/plot_tuning_curve.m
%       /media/adamX/m3ha/data_dclamp/compare_sse.m
%       /media/adamX/m3ha/data_dclamp/find_initial_slopes.m

% File History:
% 2016-12-08 Created
% 2018-06-11 Modified to use various outlier methods
% 2019-03-14 Return original data if empty
% 2019-03-14 Added 'ReplaceWithNans' as an optional argument

%% Hard-coded parameters
validOutlierMethods = {'boxplot', 'isoutlier', ...
                        'fiveStds', 'threeStds', 'twoStds'};

%% Default values for optional arguments
wl2iqrDefault = 1.5;                % same as the Matlab function boxplot()
outlierMethodDefault = 'boxplot';   % use built-in isoutlier function
replaceWithNansDefault = false;     % remove by default
plotFlagDefault = false;            % don't plot box plots by default

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
addRequired(iP, 'oldData', ... 
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'WL2IQR', wl2iqrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'OutlierMethod', outlierMethodDefault, ...
    @(x) any(validatestring(x, validOutlierMethods)));
addParameter(iP, 'ReplaceWithNans', replaceWithNansDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, oldData, varargin{:});
wl2iqr = iP.Results.WL2IQR;
outlierMethod = validatestring(iP.Results.OutlierMethod, validOutlierMethods);
replaceWithNans = iP.Results.ReplaceWithNans;
plotFlag = iP.Results.PlotFlag;

%% Preparation
% Force row vectors as column vectors
oldData = force_column_vector(oldData, 'IgnoreNonVectors', true);

%% Check if there is data
if isempty(oldData)
    newData = oldData;
    rowsToKeep = [];
    return
else
    % Initialize rowsToKeep
    nRows = size(oldData, 1);
    rowsToKeep = transpose(1:nRows);
end

%% Remove outliers
% Check whether each data point is within range
switch outlierMethod
case 'boxplot'
    % Compute the quartiles for each column
    %   Note: each row corresponds to a quartile
    Q = quantile(oldData, [0.25; 0.5; 0.75]);

    % Extract the first quartiles for each column
    q1 = Q(1, :);

    % Extract the third quartiles for each column
    q3 = Q(3, :);

    % Compute the interquartile range for each column
    IQR = q3 - q1;

    % Compute the whisker maximum
    highbar = q3 + wl2iqr * IQR;

    % Compute the whisker minimum
    lowbar = q1 - wl2iqr * IQR;

    % Check whether each data point is within the whiskers
    withinRange = oldData >= lowbar & oldData <= highbar;
case 'isoutlier'
    % Use the built-in isoutlier() function
    withinRange = ~isoutlier(oldData);
case {'fiveStds', 'threeStds', 'twoStds'}
    % Compute the mean of each column
    meanX = mean(oldData);

    % Compute the standard deviation of each column
    stdX = std(oldData);

    % Get the number of standard deviations away from the mean
    if strcmp(outlierMethod, 'fiveStds')
        nStds = 5;
    elseif strcmp(outlierMethod, 'threeStds')
        nStds = 3;
    elseif strcmp(outlierMethod, 'twoStds')
        nStds = 2;
    end

    % Compute the high threshold
    highbar = meanX + nStds * stdX;

    % Compute the low threshold
    lowbar = meanX - nStds * stdX;

    % Check whether each data point is within 
    %   nStds standard deviations of the mean 
    withinRange = oldData >= lowbar & oldData <= highbar;
end

% Remove outliers
if replaceWithNans
    % Initialize with old data
    newData = oldData;

    % Replace points not within range with NaN
    newData(~withinRange) = NaN;
else
    % Only include rows with all points within range
    rowsToKeep = find(all(withinRange, 2));

    % Extract the new data
    newData = oldData(rowsToKeep, :);
end

%% Plot boxplots for verification
if plotFlag
    figure();
    boxplot(oldData);
    figure();
    boxplot(newData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[newData, rowsToKeep] = remove_outliers (oldData, wl2iqr, plotFlag)
%% Check arguments
if nargin < 1
    error('Not enough input arguments, type ''help remove_outliers'' for usage');
elseif ~isnumeric(oldData)
    error('First argument must be a numeric array!');
elseif nargin >= 2 && (~isnumeric(wl2iqr) || length(wl2iqr) ~= 1)
    error('wl2iqr must be a single number!');
elseif nargin >= 3 && ~(plotFlag == 0 || plotFlag == 1)
    error('plotFlag must be either 0 or 1!');
end

%% Set defaults for optional arguments
if nargin < 2
    wl2iqr = 1.5;
end
if nargin < 3
    plotFlag = 0;
end

% Decide on whether each row of data should remain, initially all true
toleave = true(nRows, 1);

% Take out any rows that contain at least an outlier
for iRow = 1:nRows
    if sum([sum(oldData(iRow, :) > highbar), sum(oldData(iRow, :) < lowbar)])
        toleave(iRow) = false;
    else
        toleave(iRow) = true;
    end
end

% Find the original indices that will remain
rowsToKeep = find(toleave);

% Get the total number of data points for each group
nRows = size(oldData, 1);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
