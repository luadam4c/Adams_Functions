function testResults = test_var_difference (dataTable, measureVars, groupingVar, varargin)
%% Tests whether groups are different for each measured variable
% Usage: testResults = test_var_difference (dataTable, measureVars, groupingVar, varargin)
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
%       measureVars - name of measured variables
%                   must be a string array or a cell array of character arrays
%       groupingVar - name of categorical variable(s) to test difference on
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
%                   - Any other parameter-value pair for 
%                       the test_difference() function
%
% Requires:
%       cd/construct_fullpath.m
%       cd/convert_to_char.m
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/force_string_end.m
%       cd/match_format_vector_sets.m
%       cd/plot_histogram.m
%       cd/test_difference.m
%       cd/unique_groups.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-01-09 Created by Adam Lu
% 2019-01-10 Added 'Prefix' as an optional argument
% 2019-09-01 Renamed as test_var_difference 
% 2019-10-12 Now uses test_normality.m
% 2020-02-14 Now uses unique_groups.m
% 2020-02-14 Pulled out test_difference.m
% 

%% Hard-coded parameters
% TODO: Make these parameters
plotHistograms = true;
saveHistFlag = true; %false;    % default false
histFigNames = '';
histStyle = 'overlapped';

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'dataTable', ...
    @(x) validateattributes(x, {'table'}, {'3d'}));
addRequired(iP, 'measureVars', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['measureVars must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'groupingVar', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['groupingVar must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, dataTable, measureVars, groupingVar, varargin{:});
sheetName = iP.Results.SheetName;
prefix = iP.Results.Prefix;
outFolder = iP.Results.OutFolder;

% Keep unmatched arguments for the test_difference() function
otherArguments = iP.Unmatched;

%% Preparation
% Extract the grouping vector(s)
grouping = dataTable{:, groupingVar};

% Count the number of grouping vector(s)
nParams = size(grouping, 2);

% Force as a column cell array of character vectors
measureVars = force_column_cell(measureVars);


% Extract the y vectors
measureData = cellfun(@(x) dataTable{:, x}, ...
                        measureVars, 'UniformOutput', false);

% Creat full path to spreadsheet file if requested
if ~isempty(sheetName)
    sheetName = construct_fullpath(sheetName, 'Directory', outFolder);
end

% Make sure the prefix ends in '_'
prefix = force_string_end(prefix, '_', 'OnlyIfNonempty', true);

%% Perform the appropriate comparison test
if nParams == 1
    % Perform t test, rank sum test or ANOVA
    statsStructs = cellfun(@(a) test_difference(a, grouping, ...
                            otherArguments), measureData);

    % Convert to a table
    testResults = struct2table(statsStructs, 'AsArray', true);
else
    % Create linear models
    models = cellfun(@(a) fitlm(grouping, a), measureData, ...
                        'UniformOutput', false);

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
xVarStr = convert_to_char(groupingVar, 'SingleOutput', true, 'Delimiter', '_and_');

% Create relationship strings
relationshipStrs = strcat(measureVars, '_vs_', xVarStr);

% Use relationship strings as row names
testResults.Properties.RowNames = relationshipStrs;

%% Generate overlapped histograms
if plotHistograms && nParams == 1
    % Create figure names if not provided
    if saveHistFlag && isempty(histFigNames)
        histFigNames = strcat(prefix, 'histogram_', relationshipStrs, '.png');
        histFigNames = construct_fullpath(histFigNames, 'Directory', outFolder);
    else
        histFigNames = match_format_vector_sets(histFigNames, measureVars);
    end

    % Count the number of readout variables
    nOut = numel(measureVars);

    % Get the unique grouping values
    uniqueGroups = unique_groups(grouping, 'IgnoreNaN', true);

    % Get the unique group names
    groupNames = convert_to_char(uniqueGroups);
    groupNames = replace(groupNames, '-', 'neg');

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
    figTitles = strcat(measureVars, ": ", pValueText);

    % Create grouping labels with pNorm values attached
    groupingLabels = cellfun(@(x) strcat(groupNames, ": ", x), ...
                            pNormAvgText, 'UniformOutput', false);

    % Plot grouped histograms
    [bars, figures] = ...
        cellfun(@(x, y, z, w, v) plot_histogram(x, 'Grouping', grouping, ...
                    'GroupedStyle', histStyle, ...
                    'XLabel', y, 'GroupingLabels', z, ...
                    'FigName', w, 'FigTitle', v), ...
                    measureData, measureVars, groupingLabels, ...
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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
