function statsStruct = test_difference (data, varargin)
%% Performs the appropriate test between groups based on the normality
% Usage: statsStruct = test_difference (data, grouping (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       data1 = randn(100, 1);
%       data2 = randn(100, 1) + 1;
%       data3 = rand(100, 1);
%       data4 = rand(100, 1) + 1;
%       data5 = [data1, data2];
%       data6 = [data3, data4];
%       statsStruct1 = test_difference(data1)
%       statsStruct2 = test_difference(data2)
%       statsStruct3 = test_difference(data3)
%       statsStruct4 = test_difference(data4)
%       statsStruct5 = test_difference(data5)
%       statsStruct6 = test_difference(data6)
%
% Outputs:
%       statsStruct - a structure with fields:
%                       isDifferentXXX (between each pair of groups)
%                       testFunction
%                       isNormal
%                       pNormAvg
%                       pNormLill
%                       pNormAd
%                       pNormJb
%                   specified as a scalar structure
%
% Arguments:
%       data        - data values
%                   must be an array
%       grouping     - (opt) grouping vector
%                   must be an array
%                   default == [] (data treated as a single group)
%       varargin    - 'UniqueGroups': unique groups
%                   must be an array
%                   default == unique_groups(grouping) omitting NaNs
%                   - 'GroupNames': unique group names
%                   must be a character array, a string array
%                       or a cell array of character arrays
%                   default == converted from uniqueGroups
%                   - 'AlphaNormality': significance level for normality test
%                   must be a positive scalar
%                   default == 0.05
%                   - 'AlphaDifference': significance level for difference test
%                   must be a positive scalar
%                   default == 0.05
%                   - 'IsPaired': whether data is paired
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'DisplayAnova': whether to display ANOVA table
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/argfun.m
%       cd/compute_pairwise_differences.m
%       cd/convert_to_char.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_grouping_by_vectors.m
%       cd/test_normality.m
%       cd/unique_groups.m
%
% Used by:
%       cd/test_var_difference.m

% File History:
% 2020-02-14 Moved from test_var_difference.m
% 2020-02-14 Added the input parser
% 2020-02-14 Added 'IsPaired' as an optional argument
% 

%% Hard-coded parameters
nullMean = 0;                   % null hypothesis mean for one-sample test
nullMedian = 0;                 % null hypothesis median for one-sample test

%% Default values for optional arguments
groupingDefault  = [];          % set later
uniqueGroupsDefault = [];       % set later
groupNamesDefault = {};         % set later
alphaNormalityDefault = 0.05;   % significance level for normality test
alphaDifferenceDefault = 0.05;  % significance level for difference test
isPairedDefault = false;        % data is not paired by default
displayAnovaDefault = true;     % display ANOVA table by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'data');

% Add optional inputs to the Input Parser
addOptional(iP, 'grouping', groupingDefault);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'UniqueGroups', uniqueGroupsDefault);
addParameter(iP, 'GroupNames', groupNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['GroupNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'AlphaNormality', alphaNormalityDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'AlphaDifference', alphaDifferenceDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'IsPaired', isPairedDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'DisplayAnova', displayAnovaDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, data, varargin{:});
grouping = iP.Results.grouping;
uniqueGroups = iP.Results.UniqueGroups;
groupNames = iP.Results.GroupNames;
alphaNormality = iP.Results.AlphaNormality;
alphaDifference = iP.Results.AlphaDifference;
isPaired = iP.Results.IsPaired;
displayAnova = iP.Results.DisplayAnova;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Preparation
% Create a grouping array if not provided
if isempty(grouping)
    grouping = create_grouping_by_vectors(data);
end

% Linearize data and grouping array
data = data(:);
grouping = grouping(:);

% Get the unique grouping values
if isempty(uniqueGroups)
    uniqueGroups = unique_groups(grouping, 'IgnoreNaN', true);
end

% Get the unique group names, but 
if isempty(groupNames)
    groupNames = convert_to_char(uniqueGroups, 'ForceCellOutput', true);
end

% Replace '-' with 'neg'
%   Note: Field names of structures cannot contain '-'
groupNames = replace(groupNames, '-', 'neg');

% Count the number of groups
nGroups = numel(uniqueGroups);

% Create pNormStrs
[pNormAvgStrs, pNormLillStrs, pNormAdStrs, pNormJbStrs] = ...
    argfun(@(x) strcat(x, groupNames), ...
            'pNormAvg_', 'pNormLill_', 'pNormAd_', 'pNormJb_');

%% Do the job
% If there are too many NaNs, return
if sum(isnan(data)) >= numel(data) / 2
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
if iscell(uniqueGroups)
    data = cellfun(@(w) data(ismatch(grouping, w)), uniqueGroups, ...
                    'UniformOutput', false);
else
    data = arrayfun(@(w) data(ismatch(grouping, w)), uniqueGroups, ...
                    'UniformOutput', false);
end

% Compute pairwise differences
if isPaired
    diffData = compute_pairwise_differences(data);
end

% Test the normality of the data
if isPaired
    % Test the normality of pairwise differences
    % TODO    
else
    % Test the normality of each column
    [isNormal, pTable] = test_normality(data, 'SigLevel', alphaNormality);
    pNormLill = pTable.pNormLill; 
    pNormAd = pTable.pNormAd; 
    pNormJb = pTable.pNormJb;
    pNormAvg = pTable.pNormAvg;
end

%% Perform the correct difference test among groups
if nGroups == 1
    if all(isNormal)
        % Perform a 1-sample t-test 
        %   (tests difference of mean with nullMean)
        [isDifferent, pValue] = ttest(data{1}, nullMean, ...
                                        'Alpha', alphaDifference);

        % Store test function
        testFunction = 'ttest';
    else
        % Perform a 1-sample signed-rank test 
        %   (tests difference of median with nullMedian)
        [pValue, isDifferent] = signrank(data{1}, nullMedian, ...
                                        'Alpha', alphaDifference);

        % Store test function
        testFunction = 'signrank';
    end
elseif nGroups == 2
    if all(isNormal)
        % Perform a 2-sample t-test (tests difference between means)
        [isDifferent, pValue] = ttest2(data{1}, data{2}, ...
                                        'Alpha', alphaDifference);

        % Store test function
        testFunction = 'ttest2';
    else
        % Perform a Wilcoxon rank-sum test (tests difference between medians)
        [pValue, isDifferent] = ranksum(data{1}, data{2}, ...
                                        'Alpha', alphaDifference);

        % Store test function
        testFunction = 'ranksum';
    end
else
    % Decide whether to display the ANOVA table and pairwise comparison graphs
    if displayAnova
        displayOpt = 'on';
    else
        displayOpt = 'off';
    end

    % Perform a one-way ANOVA (tests difference among means)
    [pValue, ~, stats] = anova1(data, grouping, displayOpt);

    % Decide whether the group means are different
    isDifferent = pValue < alphaDifference;

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
        isDifferentThis = pValueThis < alphaDifference;

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%