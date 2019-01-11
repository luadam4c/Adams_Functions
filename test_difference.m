function testResults = test_difference (dataTable, yVars, xVar, varargin)
%% Tests whether groups are different for each measured variable
% Usage: testResults = test_difference (dataTable, yVars, xVar, varargin)
% Explanation:
%       TODO
% Example(s):
%       data1 = [rand(60, 1); rand(60, 1) + 0.5];
%       data2 = [randn(60, 1); randn(60, 1) + 1];
%       grps1 = [2 * ones(60, 1); 1 * ones(60, 1)];
%       grps2 = [10 * ones(60, 1); -8 * ones(60, 1)];
%       tble1 = table(data1, data2, grps1, grps2);
%       testResults1 = test_difference(tble1, 'data1', 'grps1')
%       testResults2 = test_difference(tble1, {'data1', 'data2'}, 'grps2')
%       testResults3 = test_difference(tble1, {'data1', 'data2'}, {'grps1', 'grps2'})
%       data3 = [randn(60, 1); randn(60, 1) + 1; randn(60, 1) - 1];
%       data4 = [rand(60, 1); rand(60, 1) + 1; rand(60, 1) - 1];
%       grps3 = [2 * ones(60, 1); 1 * ones(60, 1); 3 * ones(60, 1)];
%       grps4 = [7 * ones(60, 1); 8 * ones(60, 1); 6 * ones(60, 1)];
%       tble3 = table(data3, data4, grps3, grps4);
%       testResults4 = test_difference(tble3, {'data3', 'data4'}, {'grps3', 'grps4'})
% Outputs:
%       testResults - a table with each row corresponding to a measured variable
%                           and columns:
%                       isDifferent - whether the groups are different
%                       pValue      - p value for each test
%                       testFunction- test function used
%                       isNormal    - whether groups are normally distributed
%                       h           - cell array of handles to histograms
%                       histFigNames- histogram file names
%                   specified as a table
% Arguments:
%       dataTable   - data table
%                   must be a table
%       yVars       - name of measured variables
%                   must be a string array or a cell array of character arrays
%       xVar        - name of categorical variable(s) to test difference on
%                   must be a string array or a cell array of character arrays
%       varargin    - 'SheetName' - spreadsheet file name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == ''
%                   - 'OutFolder': directory to place spreadsheet file
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/argfun.m
%       cd/construct_fullpath.m
%       cd/convert_to_char.m
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/match_format_vector_sets.m
%       cd/plot_grouped_histogram.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-01-09 Created by Adam Lu
% 2019-01-10 Added 'Prefix' as an optional argument
% 

%% Hard-coded parameters
% TODO: Make these parameters
plotHistograms = true;
saveHistFlag = false;
histFigNames = '';
alphaNormality = 0.05;  % significance level for normality test
alphaTest = 0.05;       % significance level for hypothesis test
toDisplay = false;      % whether to display ANOVA table

%% Default values for optional arguments
sheetNameDefault = '';
prefixDefault = '';   	% prepend nothing to file names by default
outFolderDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'dataTable', ...
    @(x) validateattributes(x, {'table'}, {'3d'}));
addRequired(iP, 'yVars', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['yVars must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'xVar', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['yVars must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, dataTable, yVars, xVar, varargin{:});
sheetName = iP.Results.SheetName;
prefix = iP.Results.Prefix;
outFolder = iP.Results.OutFolder;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Extract the x vector(s)
xData = dataTable{:, xVar};

% Count the number of x vectors
nParams = size(xData, 2);

% Force as a column cell array of character vectors
yVars = force_column_cell(yVars);

% Extract the y vectors
yData = cellfun(@(x) dataTable{:, x}, yVars, 'UniformOutput', false);

% Creat full path to spreadsheet file if requested
if ~isempty(sheetName)
    sheetName = construct_fullpath(sheetName, 'Directory', outFolder);
end

%% Generate overlapped histograms
if plotHistograms && nParams == 1
    % Create figure names if not provided
    if saveHistFlag && isempty(histFigNames)
        histFigNames = strcat(prefix, '_histogram_', yVars, '.png');
        histFigNames = construct_fullpath(histFigNames, 'Directory', outFolder);
    else
        histFigNames = match_format_vector_sets(histFigNames, yVars);
    end

    % Plot grouped histograms
    h = cellfun(@(x, y, z) plot_grouped_histogram(x, xData, 'XLabel', y, ...
                        'FigName', z, 'Style', 'overlapped'), ...
                        yData, yVars, histFigNames, 'UniformOutput', false);
end

%% Perform the appropriate comparison test
if nParams == 1
    % Perform t test, rank sum test or ANOVA
    statsStructs = cellfun(@(y) test_for_one_variable(xData, y, ...
                        alphaNormality, alphaTest, toDisplay), yData);

    % Convert to a table
    testResults = struct2table(statsStructs);

    % Add histograms if plotted
    if plotHistograms
        testResults = addvars(testResults, histFigNames, h);
    end
else
    % Create linear models
    models = cellfun(@(y) fitlm(xData, y), yData, 'UniformOutput', false);

    % Compute ANOVA statistics tables
    results = cellfun(@(x) anova(x), models, 'UniformOutput', false);

    % Extract statistics from each table
    sumOfSquares = cellfun(@(x) x.SumSq, results, 'UniformOutput', false);
    degreeFreedom = cellfun(@(x) x.DF, results, 'UniformOutput', false);
    meanSquares = cellfun(@(x) x.MeanSq, results, 'UniformOutput', false);
    pValue = cellfun(@(x) x.pValue, results, 'UniformOutput', false);
    fStatistic = cellfun(@(x) x.F, results, 'UniformOutput', false);

    % Decide whether groups are different
    isDifferent = cellfun(@(x) x < alphaTest, pValue, 'UniformOutput', false);

    % Place results in a table
    % TODO: Better way of organizing
    testResults = table(isDifferent, pValue, sumOfSquares, degreeFreedom, ...
                        meanSquares, fStatistic);
end

% Create row names
testResults.Properties.RowNames = strcat(yVars, '_vs_', char(strjoin(xVar, '_and_')));

%% Save the table if requested
if ~isempty(sheetName)
    writetable(testResults, sheetName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function statsStruct = test_for_one_variable(xData, yData, alphaNormality, ...
                                        alphaTest, toDisplay)
%% Performs the appropriate test based on the normality

% Get the unique x values
uniqueX = unique(xData);

% Count the number of groups
nGroups = numel(uniqueX);

% Separate the data into groups
%   Note: data will become a cell array
if iscell(uniqueX)
    data = cellfun(@(w) yData(ismatch(xData, w)), uniqueX, ...
                    'UniformOutput', false);
else
    data = arrayfun(@(w) yData(ismatch(xData, w)), uniqueX, ...
                    'UniformOutput', false);
end

% Apply the Lilliefors test for normality to each group
[~, pNormalityLill] = ...
    cellfun(@(x) lillietest(x, 'Alpha', alphaNormality), data);

% Apply the Anderson-Darling test for normality to each group
[~, pNormalityAd] = cellfun(@(x) adtest(x, 'Alpha', alphaNormality), data);

% Apply the Jarque-Bera test for normality to each group
[~, pNormalityJb] = cellfun(@(x) jbtest(x, alphaNormality), data);

% Apply the One-sample Kolmogorov-Smirnov test for normality to each group
[~, pNormalityKs] = cellfun(@(x) kstest(x, 'Alpha', alphaNormality), data);

% Place all p values for normality together in a matrix
%   Note: each row is a group; each column is a different test
pNormalityMat = [pNormalityLill, pNormalityAd, pNormalityJb, pNormalityKs];

% Take the geometric mean of the p values from different tests
pNormality = compute_weighted_average(pNormalityMat, 'DimToOperate', 2, ...
                                        'AverageMethod', 'geometric');

% Normality is satified if all p values are greater than the significance level
isNormal = all(pNormality >= alphaNormality);

% Perform the correct test
if nGroups == 1
    % Perform a 1-sample t-test (tests difference of mean with 0)
    [isDifferent, pValue] = ttest1(data{:}, 'Alpha', alphaTest);

    % Store test function
    testFunction = 'ttest1';
elseif nGroups == 2 && isNormal
    % Perform a 2-sample t-test (tests difference between means)
    [isDifferent, pValue] = ttest2(data{:}, 'Alpha', alphaTest);

    % Store test function
    testFunction = 'ttest2';
elseif nGroups == 2 && ~isNormal
    % Perform a Wilcoxon rank-sum test (tests difference between medians)
    [pValue, isDifferent] = ranksum(data{:}, 'Alpha', alphaTest);

    % Store test function
    testFunction = 'ranksum';
else
    % Decide whether to display the ANOVA table
    if toDisplay
        displayOpt = 'on';
    else
        displayOpt = 'off';
    end

    % Perform a one-way ANOVA (tests difference among means)
    [pValue, ~, stats] = anova1(yData, xData, displayOpt);

    % Decide whether the group means are different
    isDifferent = pValue < alphaTest;

    % Store test function
    testFunction = 'anova1';
end

% Store in statsStruct
statsStruct.isDifferent = isDifferent;
statsStruct.pValue = pValue;
statsStruct.testFunction = testFunction;
statsStruct.isNormal = isNormal;
statsStruct.pNormalityLill = pNormalityLill;
statsStruct.pNormalityAd = pNormalityAd;
statsStruct.pNormalityJb = pNormalityJb;
statsStruct.pNormalityKs = pNormalityKs;

% Decide whether there is a difference between each pair of groups
if nGroups > 2
    % Apply multcompare()
    %   Note: Each row is a pair. Columns:
    %           1 - index of first group
    %           2 - index of second group
    %           3 - lower confidence interval of difference
    %           4 - mean difference
    %           5 - upper confidence interval of difference
    %           6 - p value
    otherStats = multcompare(stats);

    % Create a vector of first groups
    [firstGroupNames, secondGroupNames] = ...
        argfun(@(x) convert_to_char(uniqueX(otherStats(:, x))), 1, 2);

    % Create pValueStrs
    pValueStrs = strcat('pValue_', firstGroupNames, '_', secondGroupNames);

    % Count the number of pairs
    nPairs = size(otherStats, 1);
    
    % Store p values in statsStruct
    for iPair = 1:nPairs
        statsStruct.(pValueStrs{iPair}) = otherStats(iPair, 6);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if isstring(yVars)
    % Use arrayfun
    yData = arrayfun(@(x) dataTable{:, x}, yVars, 'UniformOutput', false);
else
end

%% If not compiled, add directories to search path for required functions
%       /home/Matlab/Downloaded_Functions/swtest.m
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for swtest.m
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions'));
end
% Apply the Shapiro-Wilk Test for normality
[isNotNormal, pNormality] = cellfun(@(y) swtest(y, alphaNormality), yData);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%