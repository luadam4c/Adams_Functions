function testResults = test_var_difference (dataTable, yVars, xVar, varargin)
%% Tests whether groups are different for each measured variable
% Usage: testResults = test_var_difference (dataTable, yVars, xVar, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       n = 180
%       data1 = [rand(n, 1); rand(n, 1) + 0.5];
%       data2 = [randn(n, 1); randn(n, 1) + 1];
%       data3 = [randn(n, 1); randn(n, 1)];
%       grps1 = [2 * ones(n, 1); 1 * ones(n, 1)];
%       grps2 = [10 * ones(n, 1); -8 * ones(n, 1)];
%       tble1 = table(data1, data2, data3, grps1, grps2);
%       testResults1 = test_var_difference(tble1, 'data1', 'grps1')
%       testResults2 = test_var_difference(tble1, {'data1', 'data2', 'data3'}, 'grps2')
%       testResults3 = test_var_difference(tble1, {'data1', 'data2'}, {'grps1', 'grps2'})
%       data4 = [randn(n, 1); randn(n, 1); randn(n, 1)];
%       data5 = [randn(n, 1); randn(n, 1) + 1; randn(n, 1) - 1];
%       grps4 = [2 * ones(n, 1); -1 * ones(n, 1); 3 * ones(n, 1)];
%       grps5 = [7 * ones(n, 1); 8 * ones(n, 1); 6 * ones(n, 1)];
%       tble4 = table(data4, data5, grps4, grps5);
%       testResults4 = test_var_difference(tble4, {'data4', 'data5'}, 'grps4')
%       testResults5 = test_var_difference(tble4, {'data4', 'data5'}, 'grps5')
%       testResults6 = test_var_difference(tble4, {'data4', 'data5'}, {'grps4', 'grps5'})
%
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
%
% Arguments:
%       dataTable   - data table
%                   must be a table
%       yVars       - name of measured variables
%                   must be a string array or a cell array of character arrays
%       xVar        - name of categorical variable(s) to test difference on
%                   must be a string array or a cell array of character arrays
%       varargin    - 'SheetName' - spreadsheet file name for saving
%                   must be a string scalar or a character vector
%                   default == '' (don't save)
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
%       cd/force_string_end.m
%       cd/match_format_vector_sets.m
%       cd/plot_histogram.m
%       cd/struct2arglist.m
%       cd/test_normality.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-01-09 Created by Adam Lu
% 2019-01-10 Added 'Prefix' as an optional argument
% 2019-09-01 Renamed as test_var_difference 
% 2019-10-12 Now uses test_normality.m
% TODO: Pull out test_var_difference_helper.m that will take two vectors as required arguments
% 

%% Hard-coded parameters
% TODO: Make these parameters
plotHistograms = true;
saveHistFlag = true; %false;    % default false
toDisplay = true;       % whether to display ANOVA table
histFigNames = '';
histStyle = 'overlapped';
alphaNormality = 0.05;  % significance level for normality test
alphaTest = 0.05;       % significance level for hypothesis test

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

% Get the unique x values
uniqueX = unique(xData);

% Get the unique group names as a cell array of character arrays
groupNames = convert_to_char(uniqueX);
groupNames = replace(groupNames, '-', 'neg');
if ischar(groupNames)
    groupNames = {groupNames};
end

% Force as a column cell array of character vectors
yVars = force_column_cell(yVars);

% Count the number of readout variables
nOut = numel(yVars);

% Extract the y vectors
yData = cellfun(@(x) dataTable{:, x}, yVars, 'UniformOutput', false);

% Creat full path to spreadsheet file if requested
if ~isempty(sheetName)
    sheetName = construct_fullpath(sheetName, 'Directory', outFolder);
end

% Make sure the prefix ends in '_'
prefix = force_string_end(prefix, '_', 'OnlyIfNonempty', true);

%% Perform the appropriate comparison test
if nParams == 1
    % Perform t test, rank sum test or ANOVA
    statsStructs = cellfun(@(y) test_var_difference_helper(xData, y, ...
                        uniqueX, groupNames, alphaNormality, ...
                        alphaTest, toDisplay), yData);

    % Convert to a table
    testResults = struct2table(statsStructs, 'AsArray', true);
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

% Combine all x variable strings into one string
xVarStr = convert_to_char(xVar, 'SingleOutput', true, 'Delimiter', '_and_');

% Create relationship strings
relationshipStrs = strcat(yVars, '_vs_', xVarStr);

% Use relationship strings as row names
testResults.Properties.RowNames = relationshipStrs;

%% Generate overlapped histograms
if plotHistograms && nParams == 1
    % Create figure names if not provided
    if saveHistFlag && isempty(histFigNames)
        histFigNames = strcat(prefix, 'histogram_', relationshipStrs, '.png');
        histFigNames = construct_fullpath(histFigNames, 'Directory', outFolder);
    else
        histFigNames = match_format_vector_sets(histFigNames, yVars);
    end

    % Create pNormAvgStrs
    pNormAvgStrs = strcat('pNormAvg_', groupNames);

    % Extract information for labelling
    isDifferent = testResults.isDifferent;
    pValue = testResults.pValue;
    pNormAvg = testResults{:, pNormAvgStrs};
    
    % Extract information for labelling
    pValueText = strcat("pValue = ", num2str(pValue, 1));
    pNormAvgText = arrayfun(@(x) strcat("pNorm = ", ...
                            num2str(transpose(pNormAvg(x, :)), 1)), ...
                            transpose(1:nOut), 'UniformOutput', false);
    
    % Create titles with pValueText attached
    figTitles = strcat(yVars, ": ", pValueText);

    % Create grouping labels with pNorm values attached
    groupingLabels = cellfun(@(x) strcat(groupNames, ": ", x), ...
                            pNormAvgText, 'UniformOutput', false);

    % Plot grouped histograms
    [bars, figures] = ...
        cellfun(@(x, y, z, w, v) plot_histogram(x, 'Grouping', xData, ...
                    'GroupedStyle', histStyle, ...
                    'XLabel', y, 'GroupingLabels', z, ...
                    'FigName', w, 'FigTitle', v), ...
                    yData, yVars, groupingLabels, ...
                    histFigNames, figTitles, ...
                    'UniformOutput', false);  

    % Append histogram to results
    testResults = addvars(testResults, histFigNames, bars, figures);
end

%% Save the table if requested
if ~isempty(sheetName)
    writetable(testResults, sheetName, 'WriteRowNames', true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function statsStruct = test_var_difference_helper (xData, yData, uniqueX, ...
                                        groupNames, alphaNormality, ...
                                        alphaTest, toDisplay)
%% Performs the appropriate test based on the normality

% Get the unique x values
if isempty(uniqueX)
    uniqueX = unique(xData);
end

% Get the unique group names, but 
if isempty(groupNames)
    groupNames = convert_to_char(uniqueX);
end

% Replace '-' with 'neg'
%   Note: Field names of structures cannot contain '-'
groupNames = replace(groupNames, '-', 'neg');

% Count the number of groups
nGroups = numel(uniqueX);

% Create pNormStrs
[pNormAvgStrs, pNormLillStrs, pNormAdStrs, pNormJbStrs] = ...
    argfun(@(x) strcat(x, groupNames), ...
            'pNormAvg_', 'pNormLill_', 'pNormAd_', 'pNormJb_');

% If there are too many NaNs, return
if sum(isnan(yData)) >= numel(yData) / 2
    statsStruct.isDifferent = false;
    statsStruct.pValue = NaN;
    statsStruct.testFunction = 'none';
    statsStruct.isNormal = false;
    statsStruct.pNormAvg = NaN;
    statsStruct.pNormLill = NaN;
    statsStruct.pNormAd = NaN;
    statsStruct.pNormJb = NaN;
    for iGroup = 1:nGroups
        statsStruct.(pNormAvgStrs{iGroup}) = NaN;
        statsStruct.(pNormLillStrs{iGroup}) = NaN;
        statsStruct.(pNormAdStrs{iGroup}) = NaN;
        statsStruct.(pNormJbStrs{iGroup}) = NaN;
    end
    return
end

% Separate the data into groups
%   Note: data will become a cell array
if iscell(uniqueX)
    data = cellfun(@(w) yData(ismatch(xData, w)), uniqueX, ...
                    'UniformOutput', false);
else
    data = arrayfun(@(w) yData(ismatch(xData, w)), uniqueX, ...
                    'UniformOutput', false);
end

% TODO: use test_normality.m
[isNormal, pTable] = test_normality(data, 'SigLevel', alphaNormality);
pNormLill = pTable.pNormLill; 
pNormAd = pTable.pNormAd; 
pNormJb = pTable.pNormJb;
pNormAvg = pTable.pNormAvg;

% Perform the correct test
if nGroups == 1
    % Perform a 1-sample t-test (tests difference of mean with 0)
    [isDifferent, pValue] = ttest1(data{:}, 'Alpha', alphaTest);

    % Store test function
    testFunction = 'ttest1';
elseif nGroups == 2 && all(isNormal)
    % Perform a 2-sample t-test (tests difference between means)
    [isDifferent, pValue] = ttest2(data{:}, 'Alpha', alphaTest);

    % Store test function
    testFunction = 'ttest2';
elseif nGroups == 2 && ~all(isNormal)
    % Perform a Wilcoxon rank-sum test (tests difference between medians)
    [pValue, isDifferent] = ranksum(data{:}, 'Alpha', alphaTest);

    % Store test function
    testFunction = 'ranksum';
else
    % Decide whether to display the ANOVA table and pairwise comparison graphs
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
    otherStats = multcompare(stats, 'Display', displayOpt);

    % Create group name vectors
    [firstGroupNames, secondGroupNames] = ...
        argfun(@(x) groupNames(otherStats(:, x)), 1, 2);

    % Create isDifferentStrs
    isDifferentStrs = strcat('isDifferent_', firstGroupNames, ...
                                '_', secondGroupNames);

    % Create pValueStrs
    pValueStrs = strcat('pValue_', firstGroupNames, '_', secondGroupNames);

    % Count the number of pairs
    nPairs = size(otherStats, 1);
    
    % Store p values in statsStruct
    for iPair = 1:nPairs
        % Extract pvalue
        pValueThis = otherStats(iPair, 6);

        % Test whether there is a difference between this pair
        isDifferentThis = pValueThis < alphaTest;

        % Store isDifferent for this pair
        statsStruct.(isDifferentStrs{iPair}) = isDifferentThis;

        % Store p values
        statsStruct.(pValueStrs{iPair}) = pValueThis;
    end
end

% Store other information
statsStruct.testFunction = testFunction;
statsStruct.isNormal = isNormal;
statsStruct.pNormAvg = pNormAvg;
statsStruct.pNormLill = pNormLill;
statsStruct.pNormAd = pNormAd;
statsStruct.pNormJb = pNormJb;

% Store normality test p values in statsStruct
for iGroup = 1:nGroups
    statsStruct.(pNormAvgStrs{iGroup}) = pNormAvg(iGroup);
    statsStruct.(pNormLillStrs{iGroup}) = pNormLill(iGroup);
    statsStruct.(pNormAdStrs{iGroup}) = pNormAd(iGroup);
    statsStruct.(pNormJbStrs{iGroup}) = pNormJb(iGroup);
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
[isNotNormal, pNorm] = cellfun(@(y) swtest(y, alphaNormality), yData);

% Apply the One-sample Kolmogorov-Smirnov test for normality to each group
[~, pNormKs] = cellfun(@(x) kstest(x, 'Alpha', alphaNormality), data);
statsStruct.pNormKs = pNormKs;
pNormMat = [pNormLill, pNormAd, pNormJb, pNormKs];

% Apply the Lilliefors test for normality to each group
[~, pNormLill] = ...
    cellfun(@(x) lillietest(x, 'Alpha', alphaNormality), data);

% Apply the Anderson-Darling test for normality to each group
[~, pNormAd] = cellfun(@(x) adtest(x, 'Alpha', alphaNormality), data);

% Apply the Jarque-Bera test for normality to each group
[~, pNormJb] = cellfun(@(x) jbtest(x, alphaNormality), data);

% Place all p values for normality together in a matrix
%   Note: each row is a group; each column is a different test
pNormMat = [pNormLill, pNormAd, pNormJb];

% Take the geometric mean of the p values from different tests
pNormAvg = compute_weighted_average(pNormMat, 'DimToOperate', 2, ...
                                        'AverageMethod', 'geometric');

% Normality is satified if all p values are greater than the significance level
isNormal = all(pNormAvg >= alphaNormality);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
