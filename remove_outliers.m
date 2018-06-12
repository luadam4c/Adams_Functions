function [newData, origInd] = remove_outliers (oldData, varargin)
%% Removes outliers from a data matrix and return a new matrix
% Usage: [newData, origInd] = remove_outliers (oldData, varargin)
% Example(s):
%       TODO
% Outputs:
%       newData     - data matrix with outlying data points removed
%                   specified as a numeric array
%       origInd     - row indices of original data matrix that 
%                       were left in place
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
%                   - 'PlotFlag': whether to plot box plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Used by:    
%       /home/Matlab/Adams_Functions/histogram_include_outofrange.m
%       /media/adamX/m3ha/data_dclamp/compare_sse.m
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m
%
% 2016-12-08 Created
% 2018-06-11 Modified to use various outlier methods
% 

%% Hard-coded parameters
validOutlierMethods = {'boxplot', 'isoutlier', ...
                        'fiveStds', 'threeStds', 'twoStds'};

%% Default values for optional arguments
wl2iqrDefault = 1.5;                % same as the Matlab function boxplot()
outlierMethodDefault = 'boxplot';   % use built-in isoutlier function
plotFlagDefault = false;            % whether to plot box plots by default

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
                                'datetime', 'duration'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'WL2IQR', wl2iqrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'OutlierMethod', outlierMethodDefault, ...
    @(x) any(validatestring(x, validOutlierMethods)));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, oldData, varargin{:});
wl2iqr = iP.Results.WL2IQR;
outlierMethod = validatestring(iP.Results.OutlierMethod, validOutlierMethods);
plotFlag = iP.Results.PlotFlag;

%% Remove outliers
switch outlierMethod
case 'boxplot'
    % Get the total number of data points
    ndp = size(oldData, 1);

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
    highbar = q3 + wl2iqr * (q3 - q1);

    % Compute the whisker minimum
    lowbar = q1 - wl2iqr * (q3 - q1);

    % Decide on whether to leave data point or not, initially all true
    toleave = ones(1, ndp);
    for k = 1:ndp
        if sum([sum(oldData(k, :) > highbar), sum(oldData(k, :) < lowbar)])
            toleave(k) = false;
        else
            toleave(k) = true;
        end
    end
    origInd = find(toleave);
    newData = oldData(origInd, :);
case 'isoutlier'
    % Take out values with the built-in isoutlier() function
    origInd = find(~isoutlier(oldData));
    newData = oldData(origInd);
case {'fiveStds', 'threeStds', 'twoStds'}
    % Compute the mean
    meanX = mean(oldData);

    % Compute the standard deviation
    stdX = std(oldData);

    % Get the number of standard deviations away from the mean
    if strcmp(outlierMethod, 'fiveStds')
        nStds = 5;
    elseif strcmp(outlierMethod, 'threeStds')
        nStds = 3;
    elseif strcmp(outlierMethod, 'twoStds')
        nStds = 2;
    end

    % Only include the points within a nStds standard deviations of the mean
    origInd = find(oldData >= meanX - nStds * stdX & ...
                    oldData <= meanX + nStds * stdX);
    newData = oldData(origInd);
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

[newData, origInd] = remove_outliers (oldData, wl2iqr, plotFlag)
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

%}