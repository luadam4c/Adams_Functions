function statsStruct = test_difference (data, varargin)
%% Performs the appropriate test between groups based on the normality
% Usage: statsStruct = test_difference (data, grouping (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       n1 = 100; n2 = 100;
%       data1 = randn(n1, 1);
%       data2 = randn(n1, 1) + 1;
%       data3 = randn(n1, 1);
%       data4 = randn(n1, 1);
%       data6 = random('Exponential', 1, n2, 1) - 1;
%       data7 = random('Exponential', 1, n2, 1) + 1;
%       data8 = random('Exponential', 1, n2, 1) - 1;
%       data9 = random('Exponential', 1, n2, 1) - 1;
%       data11 = [data1, data2];
%       data12 = [data1, data3];
%       data13 = [data6, data7];
%       data14 = [data6, data8];
%       data15 = [data1, data1 + data3];
%       data16 = [data1, data2, data3];
%       data17 = [data1, data3, data4];
%       data18 = [data6, data7, data8];
%       data19 = [data6, data8, data9];
%       data20 = [(1:n1)', (1:n1)' + data2, (1:n1)' + data3];
%       data21 = [(1:n1)', (1:n1)' + data3, (1:n1)' + data4];
%       data22 = [(1:n2)', (1:n2)' + data7, (1:n2)' + data8];
%       data23 = [(1:n2)', (1:n2)' + data8, (1:n2)' + data9];
%       statsStruct1 = test_difference(data1)
%       statsStruct2 = test_difference(data2)
%       statsStruct6 = test_difference(data6)
%       statsStruct7 = test_difference(data7)
%       statsStruct111 = test_difference(data11)
%       statsStruct112 = test_difference(data11, 'IsPaired', true)
%       statsStruct121 = test_difference(data12)
%       statsStruct122 = test_difference(data12, 'IsPaired', true)
%       statsStruct131 = test_difference(data13)
%       statsStruct132 = test_difference(data13, 'IsPaired', true)
%       statsStruct141 = test_difference(data14)
%       statsStruct142 = test_difference(data14, 'IsPaired', true)
%       statsStruct151 = test_difference(data15)
%       statsStruct152 = test_difference(data15, 'IsPaired', true)
%       statsStruct161 = test_difference(data16)
%       statsStruct162 = test_difference(data16, 'IsPaired', true)
%       statsStruct171 = test_difference(data17)
%       statsStruct172 = test_difference(data17, 'IsPaired', true)
%       statsStruct181 = test_difference(data18)
%       statsStruct182 = test_difference(data18, 'IsPaired', true)
%       statsStruct191 = test_difference(data19)
%       statsStruct192 = test_difference(data19, 'IsPaired', true)
%       statsStruct201 = test_difference(data20)
%       statsStruct202 = test_difference(data20, 'IsPaired', true)
%       statsStruct211 = test_difference(data21)
%       statsStruct212 = test_difference(data21, 'IsPaired', true)
%       statsStruct221 = test_difference(data22)
%       statsStruct222 = test_difference(data22, 'IsPaired', true)
%       statsStruct231 = test_difference(data23)
%       statsStruct232 = test_difference(data23, 'IsPaired', true)
%
%
% Outputs:
%       statsStruct - a structure with fields:
%                       symbol
%                       isDifferent
%                       pValue
%                       nSamples
%                       degreesOfFreedom
%                       testFunction
%                       symbol_Group1_Group2, etc.
%                       isDifferent_Group1_Group2, etc.
%                       pValue_Group1_Group2, etc.
%                       diffValue_Group1_Group2, etc.
%                       isNormal_Group1, etc.
%                       pNormAvg_Group1, etc.
%                       pNormLill_Group1, etc.
%                       pNormAd_Group1, etc.
%                       pNormJb_Group1, etc.
%                   specified as a scalar structure
%
% Arguments:
%       data        - data values
%                   must be an array
%       grouping     - (opt) corresponding grouping values
%                   must be an array
%                   default == create_grouping_by_vectors(data);
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
%                   - 'AlphaTwoStars': significance level for having two stars
%                                       as the symbol
%                   must be a positive scalar
%                   default == 0.01
%                   - 'AlphaThreeStars': significance level for 
%                                           having three stars as the symbol
%                   must be a positive scalar
%                   default == 0.001
%                   - 'IsPaired': whether data is paired
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'DisplayAnova': whether to display ANOVA table
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveFlag': whether to save results in a text file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'FileBase': file base name for saving
%                   must be a string scalar or a character vector
%                   default == combine_strings(groupNames)
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/argfun.m
%       cd/combine_strings.m
%       cd/compute_stats.m
%       cd/convert_to_char.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/create_grouping_by_vectors.m
%       cd/create_labels_from_numbers.m
%       cd/force_column_cell.m
%       cd/force_data_as_matrix.m
%       cd/print_structure.m
%       cd/test_normality.m
%       cd/unique_groups.m
%
% Used by:
%       cd/compare_events_pre_post_stim.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_simulate_population.m
%       cd/test_var_difference.m
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt\virt_moore.m

% File History:
% 2020-02-14 Moved from test_var_difference.m
% 2020-02-14 Added the input parser
% 2020-02-14 Added 'IsPaired' as an optional argument
% 2020-02-15 Finished implementing 'IsPaired'
% 2020-02-15 Added nonparametric tests for multiple groups
% 2020-02-16 Added meanDiff
% 2020-02-23 Now uses force_data_as_matrix.m
% 2020-02-23 Added percChange
% 2020-02-23 diffValue is now the negative of before
% 2020-02-23 Added stderrDiff, lowerDiff, upperDiff
% 2020-03-10 Added symbol, alphaTwoStars, alphaThreeStars
% 2020-05-14 Added 'SaveFlag' as an optional argument
% 2020-05-14 Added 'FileBase' as an optional argument
% TODO: 2020-03-10 Added alternative tests

%% Hard-coded parameters
nullMean = 0;                   % null hypothesis mean for one-sample test
nullMedian = 0;                 % null hypothesis median for one-sample test
diffDelimiter = '_minus_';
ranovaTimeVar = 'Time';
ranovaRowOfInterest = ['(Intercept):', ranovaTimeVar];
outSuffix = 'stats';
outExtension = 'txt';

%% Default values for optional arguments
groupingDefault  = [];          % set later
uniqueGroupsDefault = [];       % set later
groupNamesDefault = {};         % set later
alphaNormalityDefault = 0.05;   % significance level for normality test
alphaDifferenceDefault = 0.05;  % significance level for difference test
alphaTwoStarsDefault = 0.01;    % significance level for two stars
alphaThreeStarsDefault = 0.001; % significance level for three stars
isPairedDefault = false;        % data is not paired by default
displayAnovaDefault = true;     % display ANOVA table by default
saveFlagDefault = false;        % don't save results by default
fileBaseDefault = '';           % set later

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
addParameter(iP, 'AlphaTwoStars', alphaTwoStarsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'AlphaThreeStars', alphaThreeStarsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'IsPaired', isPairedDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'DisplayAnova', displayAnovaDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, data, varargin{:});
grouping = iP.Results.grouping;
uniqueGroups = iP.Results.UniqueGroups;
groupNames = iP.Results.GroupNames;
alphaNormality = iP.Results.AlphaNormality;
alphaDifference = iP.Results.AlphaDifference;
alphaTwoStars = iP.Results.AlphaTwoStars;
alphaThreeStars = iP.Results.AlphaThreeStars;
isPaired = iP.Results.IsPaired;
displayAnova = iP.Results.DisplayAnova;
saveFlag = iP.Results.SaveFlag;
fileBase = iP.Results.FileBase;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Preparation
% Force the data into a numeric matrix
dataMatrixOrig = force_data_as_matrix(data);
groupingMatrixOrig = force_data_as_matrix(grouping);

% Create a grouping matrix if not provided
if isempty(groupingMatrixOrig)
    groupingMatrixOrig = create_grouping_by_vectors(dataMatrixOrig);
end

% Remove any row with NaN if is paired
if isPaired
    rowsToKeep = ~any(isnan(dataMatrixOrig), 2);
    dataMatrix = dataMatrixOrig(rowsToKeep, :);
    groupingMatrix = groupingMatrixOrig(rowsToKeep, :);
else
    dataMatrix = dataMatrixOrig;
    groupingMatrix = groupingMatrixOrig;
end

% If groupingMatrix is empty, no data remains
%   Use original grouping matrix to generate default fields
if isempty(groupingMatrix)
    noDataRemains = true;
    groupingMatrix = groupingMatrixOrig;
else
    noDataRemains = false;
end

% Linearize data and grouping array
dataVec = dataMatrix(:);
groupingVec = groupingMatrix(:);

% Remove NaN values if not paired
if ~isPaired
    indToKeep = ~isnan(dataVec);
    dataVec = dataVec(indToKeep);
    groupingVec = groupingVec(indToKeep);
end

% Get the unique grouping values
if isempty(uniqueGroups)
    uniqueGroups = unique_groups(groupingVec, 'IgnoreNaN', true);
end

% Get the unique group names
if isempty(groupNames)
    if isnumeric(uniqueGroups)
        groupNames = create_labels_from_numbers(uniqueGroups, ...
                                                'Prefix', 'Group');
    else
        groupNames = convert_to_char(uniqueGroups, 'ForceCellOutput', true);
    end
else
    groupNames = force_column_cell(groupNames);
end

% Replace '-' with 'neg'
%   Note: Field names of structures cannot contain '-'
groupNames = replace(groupNames, '-', 'neg');

% Count the number of groups
nGroups = numel(uniqueGroups);

% Count the number of pairs
if nGroups > 1
    nPairs = nchoosek(nGroups, 2);
end

% Create group names for normality tests
if isPaired
    if nGroups > 2
        normGroupNames = strcat(groupNames, '_minus_sampleMean');
    else
        normGroupNames = {[groupNames{2}, diffDelimiter, groupNames{1}]};
    end
else
    normGroupNames = groupNames;
end

% Count the number of groups for normality tests
nNormGroups = numel(normGroupNames);

% Create strings for saving
[isNormalStrs, pNormAvgStrs, pNormLillStrs, pNormAdStrs, pNormJbStrs] = ...
    argfun(@(x) strcat(x, normGroupNames), ...
            'isNormal_', 'pNormAvg_', 'pNormLill_', 'pNormAd_', 'pNormJb_');

% Create group name vectors
if nGroups > 2
    groupNamesPaired = nchoosek(groupNames, 2);
    firstGroupNames = groupNamesPaired(:, 1);
    secondGroupNames = groupNamesPaired(:, 2);

    % Create strings
    [isDifferentStrs, symbolStrs, pValueStrs, ...
            diffValueStrs, stderrValueStrs, ...
            errValueStrs, lowerValueStrs, upperValueStrs, ...
            meanDiffStrs, stderrDiffStrs, ...
            errDiffStrs, lowerDiffStrs, upperDiffStrs, ...
            percChangeStrs] = ...
        argfun(@(a) strcat(a, firstGroupNames, '_', secondGroupNames), ...
                'isDifferent_', 'symbol_', 'pValue_', ...
                'diffValue_', 'stderrValue_', ...
                'errValue_', 'lowerValue_', 'upperValue_', ...
                'meanDiff_', 'stderrDiff_', ...
                'errDiff_', 'lowerDiff_', 'upperDiff_', ...
                'percChange_');
end

% Decide whether to display the ANOVA table and pairwise comparison graphs
if displayAnova
    displayOpt = 'on';
else
    displayOpt = 'off';
end

%% Do the job
% Initialize statsStruct
statsStruct.isDifferent = false;
statsStruct.symbol = '';
statsStruct.pValue = NaN;
statsStruct.meanDiff = NaN;
statsStruct.stderrDiff = NaN;
statsStruct.errDiff = NaN;
statsStruct.lowerDiff = NaN;
statsStruct.upperDiff = NaN;
statsStruct.percChange = NaN;
statsStruct.nSamples = NaN;
statsStruct.degreesOfFreedom = NaN;
statsStruct.testFunction = 'none';
if nGroups > 2
    for iPair = 1:nPairs
        statsStruct.(isDifferentStrs{iPair}) = NaN;
        statsStruct.(symbolStrs{iPair}) = '';
        statsStruct.(pValueStrs{iPair}) = NaN;
        statsStruct.(diffValueStrs{iPair}) = NaN;
        statsStruct.(stderrValueStrs{iPair}) = NaN;
        statsStruct.(errValueStrs{iPair}) = NaN;
        statsStruct.(lowerValueStrs{iPair}) = NaN;
        statsStruct.(upperValueStrs{iPair}) = NaN;
        statsStruct.(meanDiffStrs{iPair}) = NaN;
        statsStruct.(stderrDiffStrs{iPair}) = NaN;
        statsStruct.(errDiffStrs{iPair}) = NaN;
        statsStruct.(lowerDiffStrs{iPair}) = NaN;
        statsStruct.(upperDiffStrs{iPair}) = NaN;
        statsStruct.(percChangeStrs{iPair}) = NaN;
    end
end
for iGroup = 1:nNormGroups
    statsStruct.(isNormalStrs{iGroup}) = NaN;
    statsStruct.(pNormAvgStrs{iGroup}) = NaN;
    statsStruct.(pNormLillStrs{iGroup}) = NaN;
    statsStruct.(pNormAdStrs{iGroup}) = NaN;
    statsStruct.(pNormJbStrs{iGroup}) = NaN;
end

% Return if no data remains or if nGroups is less than 1
if noDataRemains || nGroups < 1
    return;
end

% Separate the data into groups
%   Note: data will become a cell array
if iscell(uniqueGroups)
    dataCell = cellfun(@(w) dataVec(ismatch(groupingVec, w)), ...
                        uniqueGroups, 'UniformOutput', false);
else
    dataCell = arrayfun(@(w) dataVec(ismatch(groupingVec, w)), ...
                        uniqueGroups, 'UniformOutput', false);
end

% Compute number of samples
nSamples = min(count_samples(dataCell));
statsStruct.nSamples = nSamples;

% If there is less than 2 samples, return
if nSamples < 2
    return
end

% Decide on the data for normality tests
if isPaired
    % Eliminate between-sample differences
    if nGroups > 2
        % Compute the means for each sample
        sampleMeans = nanmean(dataMatrix, 2);

        % Use the difference to the means for normalized Data
        normData = dataMatrix - sampleMeans;
    else
        % Compute pairwise differences
        normData = dataCell{2} - dataCell{1};
    end
else
    normData = dataCell;
end

% Test the normality of the data
[isNormal, pTable] = test_normality(normData, 'SigLevel', alphaNormality);
pNormLill = pTable.pNormLill; 
pNormAd = pTable.pNormAd; 
pNormJb = pTable.pNormJb;
pNormAvg = pTable.pNormAvg;

%% Compute the mean difference
if nGroups == 1
    [statsStruct.meanDiff, statsStruct.stderrDiff, ...
            statsStruct.errDiff, statsStruct.lowerDiff, statsStruct.upperDiff] = ...
        argfun(@(a) compute_stats(dataCell{1}, a, 'IgnoreNan', true), ...
                'mean', 'stderr', 'err', 'lower95', 'upper95');
    statsStruct.percChange = NaN;
elseif nGroups == 2
    meanGroup1 = nanmean(dataCell{1});
    meanGroup2 = nanmean(dataCell{2});
    percChange = ((meanGroup2 - meanGroup1) / meanGroup1) * 100;
    if isPaired
        diffData = dataCell{2} - dataCell{1};
        [statsStruct.meanDiff, statsStruct.stderrDiff, ...
                statsStruct.errDiff, statsStruct.lowerDiff, statsStruct.upperDiff] = ...
            argfun(@(a) compute_stats(diffData, a, 'IgnoreNan', true), ...
                    'mean', 'stderr', 'err', 'lower95', 'upper95');
    else
        statsStruct.meanDiff = meanGroup2 - meanGroup1;
        statsStruct.stderrDiff = NaN;
        statsStruct.errDiff = NaN;
        statsStruct.lowerDiff = NaN;
        statsStruct.upperDiff = NaN;
    end
    statsStruct.percChange = percChange;
end

%% Perform the correct difference test among groups
if nGroups == 1
    if all(isNormal)
        % Perform a 1-sample t-test 
        %   (tests difference of mean with nullMean)
        [isDifferent, pValue, ~, testStats] = ...
            ttest(dataCell{1}, nullMean, 'Alpha', alphaDifference);

        % Store test function
        testFunction = 'ttest';
    else
        % Perform a 1-sample signed-rank test 
        %   (tests difference of median with nullMedian)
        [pValue, isDifferent] = signrank(dataCell{1}, nullMedian, ...
                                        'Alpha', alphaDifference);

        % Store test function
        testFunction = 'signrank';
    end
elseif nGroups == 2
    if all(isNormal) && isPaired
        % Perform a 2-sample paired t-test (tests difference between means)
        [isDifferent, pValue, ~, testStats] = ...
            ttest(dataCell{1}, dataCell{2}, 'Alpha', alphaDifference);

        % Store test function
        testFunction = 'ttest';
    elseif all(isNormal) && ~isPaired
        % Perform a 2-sample unpaired t-test (tests difference between means)
        [isDifferent, pValue, ~, testStats] = ...
            ttest2(dataCell{1}, dataCell{2}, 'Alpha', alphaDifference);

        % Store test function
        testFunction = 'ttest2';
    elseif ~all(isNormal) && isPaired
        % Perform a Wilcoxon signed-rank test (tests difference between medians)
        [pValue, isDifferent] = ...
            signrank(dataCell{1}, dataCell{2}, 'Alpha', alphaDifference);

        % Store test function
        testFunction = 'signrank';
    elseif ~all(isNormal) && ~isPaired
        % Perform a Wilcoxon rank-sum test (tests difference between medians)
        [pValue, isDifferent] = ...
            ranksum(dataCell{1}, dataCell{2}, 'Alpha', alphaDifference);

        % Store test function
        testFunction = 'ranksum';
    end
else
    % Use the appropriate test to compare across groups
    if all(isNormal) && isPaired
        % Create a data table
        dataTable = array2table(dataMatrix, 'VariableNames', groupNames);

        % Generate a model specification string
        %   Note: Since there is no between-subject variable, use '1'
        varList = combine_strings('SubStrings', groupNames, 'Delimiter', ',');
        modelSpec = [varList, '~1'];

        % Fit a repeated-measures linear model
        rm = fitrm(dataTable, modelSpec);

        % Compute the repeated measures analysis of variance table
        ranovatbl = ranova(rm);

        % Extract the p value for testing any difference in the means
        %   across the within-subjects factors
        pValue = ranovatbl{ranovaRowOfInterest, 'pValue'};

        % Decide whether the group means are different
        isDifferent = pValue < alphaDifference;
        
        % Store test function
        testFunction = 'ranova';
    elseif all(isNormal) && ~isPaired
        % Perform a one-way ANOVA 
        %   (tests whether means of normal distributions are the same,
        %       assuming common variance)
        [pValue, ~, testStats] = anova1(dataVec, groupingVec, displayOpt);

        % Decide whether the group means are different
        isDifferent = pValue < alphaDifference;

        % Store test function
        testFunction = 'anova1';
    elseif ~all(isNormal) && isPaired
        % If there are less than two rows left, return empty
        if size(dataMatrix, 1) < 2
            % Cannot perform Friedman's test
            isDifferent = NaN;
            pValue = NaN;
        else
            % Perform the nonparametric Friedman's test 
            %   (tests whether medians are the same,
            %       assuming common continuous distributions)
            [pValue, ~, testStats] = friedman(dataMatrix, 1, displayOpt);

            % Decide whether the group medians are different
            isDifferent = pValue < alphaDifference;
        end

        % Store test function
        testFunction = 'friedman';
    elseif ~all(isNormal) && ~isPaired
        % Perform a Kruskal-Wallis test 
        %   (tests whether medians are the same,
        %       assuming common continuous distributions)
        [pValue, ~, testStats] = ...
            kruskalwallis(dataVec, groupingVec, displayOpt);

        % Decide whether the group medians are different
        isDifferent = pValue < alphaDifference;

        % Store test function
        testFunction = 'kruskalwallis';
    end
end

% Decide on the symbol
symbol = decide_on_symbol(pValue, alphaDifference, ...
                            alphaTwoStars, alphaThreeStars);

% Extract or compute degrees of freedom
switch testFunction
    case {'ttest', 'ttest2', 'anova1'}
        degreesOfFreedom = testStats.df;
    case {'signrank', 'ranksum', 'friedman', 'kruskalwallis'}
        degreesOfFreedom = NaN;
    case 'ranova'
        degreesOfFreedom = ranovatbl{ranovaRowOfInterest, 'DF'};
    otherwise
        error('testFunction unrecognized!');
end

% Store overall differences in statsStruct
statsStruct.isDifferent = isDifferent;
statsStruct.symbol = symbol;
statsStruct.pValue = pValue;
statsStruct.degreesOfFreedom = degreesOfFreedom;
statsStruct.testFunction = testFunction;

% Decide whether there is a difference between each pair of groups
if nGroups > 2 && ~isnan(pValue)
    if strcmp(testFunction, 'ranova')
        % Multiple comparison of estimated marginal means
        %   Note: otherStats is an array with 
        %               each row corresponding to a pair and columns:
        %           Time_1      - index of first group
        %           Time_2      - index of second group
        %           Difference  - difference of means
        %           StdErr      - standard error of difference
        %           pValue      - p value
        %           Lower       - lower confidence interval of difference
        %           Upper       - upper confidence interval of difference
        otherStats = multcompare(rm, ranovaTimeVar);

        % Extract columns
        firstGroupIndicesAll = otherStats{:, [ranovaTimeVar, '_1']};
        secondGroupIndicesAll = otherStats{:, [ranovaTimeVar, '_2']};
        diffValueEachPairAll = -otherStats{:, 'Difference'};
        stderrValueEachPairAll = otherStats{:, 'StdErr'};
        lowerValueEachPairAll = -otherStats{:, 'Upper'};
        upperValueEachPairAll = -otherStats{:, 'Lower'};
        pValuesEachPairAll = otherStats{:, 'pValue'};

        % Remove duplicate rows
        rowsToKeep = firstGroupIndicesAll < secondGroupIndicesAll;
        firstGroupIndices = firstGroupIndicesAll(rowsToKeep);
        secondGroupIndices = secondGroupIndicesAll(rowsToKeep);
        diffValueEachPair = diffValueEachPairAll(rowsToKeep);
        stderrValueEachPair = stderrValueEachPairAll(rowsToKeep);
        lowerValueEachPair = lowerValueEachPairAll(rowsToKeep);
        upperValueEachPair = upperValueEachPairAll(rowsToKeep);
        pValuesEachPair = pValuesEachPairAll(rowsToKeep);

        % Check the number of pairs
        if nPairs ~= numel(pValuesEachPair)
            error('Code logic error!');
        end
    else
        % Apply multcompare()
        %   Note: otherStats is a numeric matrix with 
        %               each row corresponding to a pair and columns:
        %           1 - index of first group
        %           2 - index of second group
        %           3 - lower confidence interval of difference
        %           4 - mean difference
        %           5 - upper confidence interval of difference
        %           6 - p value
        otherStats = multcompare(testStats, 'Display', displayOpt);

        % Extract columns
        firstGroupIndices = otherStats(:, 1);
        secondGroupIndices = otherStats(:, 2);
        diffValueEachPair = -otherStats(:, 4);
        lowerValueEachPair = -otherStats(:, 5);
        upperValueEachPair = -otherStats(:, 3);
        stderrValueEachPair = (upperValueEachPair - diffValueEachPair) / ...
                                    tinv(0.975, nSamples);
        pValuesEachPair = otherStats(:, 6);

        % Check the number of pairs
        if nPairs ~= size(otherStats, 1)
            error('Code logic error!');
        end
    end

    % Compute error margins for each pair
    errValueEachPair = upperValueEachPair - diffValueEachPair;

    % Compute difference data for each pair
    if isPaired
        diffDataEachPair = ...
            arrayfun(@(a, b) dataCell{b} - dataCell{a}, ...
                    firstGroupIndices, secondGroupIndices, 'UniformOutput', false);
    end
    
    % Compute mean differences for each pair
    if isPaired                
        [meanDiffEachPair, stderrDiffEachPair, ...
                errDiffEachPair, lowerDiffEachPair, upperDiffEachPair] = ...
            argfun(@(a) cellfun(@(b) compute_stats(b, a, 'IgnoreNan', true), ...
                                diffDataEachPair), ...
                    'mean', 'stderr', 'err', 'lower95', 'upper95');
    else
        meanDiffEachPair = ...
            arrayfun(@(a, b) nanmean(dataCell{b}) - nanmean(dataCell{a}), ...
                    firstGroupIndices, secondGroupIndices);
        stderrDiffEachPair = nan(nPairs, 1);
        errDiffEachPair = nan(nPairs, 1);
        lowerDiffEachPair = nan(nPairs, 1);
        upperDiffEachPair = nan(nPairs, 1);
    end

    % Compute percentage change for each pair
    percChangeEachPair = ...
        arrayfun(@(a, b) 100 * (nanmean(dataCell{b}) - nanmean(dataCell{a})) / ...
                            nanmean(dataCell{a}), ...
                firstGroupIndices, secondGroupIndices);

    % Create group name vectors
    [firstGroupNames, secondGroupNames] = ...
        argfun(@(x) groupNames(x), firstGroupIndices, secondGroupIndices);

    % Create strings again
    [isDifferentStrs, symbolStrs, pValueStrs, ...
            diffValueStrs, stderrValueStrs, ...
            errValueStrs, lowerValueStrs, upperValueStrs, ...
            meanDiffStrs, stderrDiffStrs, ...
            errDiffStrs, lowerDiffStrs, upperDiffStrs, ...
            percChangeStrs] = ...
        argfun(@(a) strcat(a, firstGroupNames, '_', secondGroupNames), ...
                'isDifferent_', 'symbol_', 'pValue_', ...
                'diffValue_', 'stderrValue_', ...
                'errValue_', 'lowerValue_', 'upperValue_', ...
                'meanDiff_', 'stderrDiff_', ...
                'errDiff_', 'lowerDiff_', 'upperDiff_', ...
                'percChange_');

    % Store p values in statsStruct
    for iPair = 1:nPairs
        % Extract pvalue
        pValueThis = pValuesEachPair(iPair);
        diffValueThis = diffValueEachPair(iPair);
        stderrValueThis = stderrValueEachPair(iPair);
        errValueThis = errValueEachPair(iPair);
        lowerValueThis = lowerValueEachPair(iPair);
        upperValueThis = upperValueEachPair(iPair);
        meanDiffThis = meanDiffEachPair(iPair);
        stderrDiffThis = stderrDiffEachPair(iPair);
        errDiffThis = errDiffEachPair(iPair);
        lowerDiffThis = lowerDiffEachPair(iPair);
        upperDiffThis = upperDiffEachPair(iPair);
        percChangeThis = percChangeEachPair(iPair);

        % Test whether there is a difference between this pair
        isDifferentThis = pValueThis < alphaDifference;

        % Decide on a symbol between this pair
        symbolThis = decide_on_symbol(pValueThis, alphaDifference, ...
                                    alphaTwoStars, alphaThreeStars);

        % Store values for this pair
        statsStruct.(isDifferentStrs{iPair}) = isDifferentThis;
        statsStruct.(symbolStrs{iPair}) = symbolThis;
        statsStruct.(pValueStrs{iPair}) = pValueThis;
        statsStruct.(diffValueStrs{iPair}) = diffValueThis;
        statsStruct.(stderrValueStrs{iPair}) = stderrValueThis;
        statsStruct.(errValueStrs{iPair}) = errValueThis;
        statsStruct.(lowerValueStrs{iPair}) = lowerValueThis;
        statsStruct.(upperValueStrs{iPair}) = upperValueThis;
        statsStruct.(meanDiffStrs{iPair}) = meanDiffThis;
        statsStruct.(stderrDiffStrs{iPair}) = stderrDiffThis;
        statsStruct.(errDiffStrs{iPair}) = errDiffThis;
        statsStruct.(lowerDiffStrs{iPair}) = lowerDiffThis;
        statsStruct.(upperDiffStrs{iPair}) = upperDiffThis;
        statsStruct.(percChangeStrs{iPair}) = percChangeThis;
    end
end

% Store normality test results in statsStruct
for iGroup = 1:nNormGroups
    statsStruct.(isNormalStrs{iGroup}) = isNormal(iGroup);
    statsStruct.(pNormAvgStrs{iGroup}) = pNormAvg(iGroup);
    statsStruct.(pNormLillStrs{iGroup}) = pNormLill(iGroup);
    statsStruct.(pNormAdStrs{iGroup}) = pNormAd(iGroup);
    statsStruct.(pNormJbStrs{iGroup}) = pNormJb(iGroup);
end

%% Save results
% Construct output file path
if saveFlag
    % Set default file base
    if isempty(fileBase)
        fileBase = combine_strings(groupNames);
    end

    % Construct full path
    filePath = strcat(fileBase, '_', outSuffix, '.', outExtension);

    % Open file for writing
    fid = fopen(filePath, 'w');

    % Print the structure
    print_structure(statsStruct, 'FileID', fid);

    % Close the file
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function symbol = decide_on_symbol(pValue, alphaDifference, ...
                                    alphaTwoStars, alphaThreeStars)
% Decide on the difference symbol

if pValue < alphaThreeStars
    symbol = '***';
elseif pValue < alphaTwoStars
    symbol = '**';
elseif pValue < alphaDifference
    symbol = '*';
else
    symbol = 'NS';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of samples
nSamples = size(dataMatrix, 1);
% Create sample numbers
sampleNumber = transpose(1:nSamples);
% Add a variable for sample number
dataTable = addvars(dataTable, sampleNumber);
% Generate a model specification string
varList = combine_strings('SubStrings', groupNames, 'Delimiter', ',');
modelSpec = [varList, '~sampleNumber'];

% List all pairs of group names
groupNamesPaired = ...
    force_column_cell(transpose(nchoosek(groupNames, 2)), ...
                        'ToLinearize', false);
% Combine strings with diffDelimiter
normGroupNames = ...
    cellfun(@(x) [x{1}, diffDelimiter, x{2}], ...
            groupNamesPaired, 'UniformOutput', false);
% Compute pairwise differences
normData = compute_pairwise_differences(dataCell);

% Compute the repeated measures analysis of variance table
[ranovatbl, betweenSubjectsSpec, ...
    withinSubjectsSpec, hypothesisValue] = ranova(rm);

% Store vectors
statsStruct.isNormal = false;
statsStruct.pNormAvg = NaN;
statsStruct.pNormLill = NaN;
statsStruct.pNormAd = NaN;
statsStruct.pNormJb = NaN;
statsStruct.isNormal = isNormal;
statsStruct.pNormAvg = pNormAvg;
statsStruct.pNormLill = pNormLill;
statsStruct.pNormAd = pNormAd;
statsStruct.pNormJb = pNormJb;

% If there are too many NaNs, return
if sum(isnan(dataVec)) >= numel(dataVec) / 2
end

dataMatrixOrig = force_matrix(data);
groupingMatrixOrig = force_matrix(grouping);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
