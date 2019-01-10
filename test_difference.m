function testResults = test_difference (dataTable, yVars, xVar, varargin)
%% Tests whether groups are different for each measured variable
% Usage: testResults = test_difference (dataTable, yVars, xVar, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       testResults - a table with each row corresponding to a measured variable
%                           and columns:
%                       isDifferent - whether the groups are different 
%                       pValue      - p value for each test
%                       h           - cell array of handles to histograms
%                       histFigNames- histogram file names
%                   specified as a table
% Arguments:
%       dataTable   - data table
%                   must be a table
%       yVars       - name of measured variables
%                   must be a string array or a cell array of character arrays
%       xVar        - name of categorical variable to test difference on
%                   must be a string or a character array
%       varargin    - 'SheetName' - spreadsheet file name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'OutFolder': directory to place spreadsheet file
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/construct_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/plot_grouped_histogram.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-01-09 Created by Adam Lu
% 

%% Hard-coded parameters
% TODO: Make these parameters
plotHistograms = true;
saveHistFlag = true;
histFigNames = '';

%% Default values for optional arguments
sheetNameDefault = '';
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
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, dataTable, varargin{:});
sheetName = iP.Results.SheetName;
outFolder = iP.Results.OutFolder;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Extract the x vector
xData = dataTable{:, xVar};

% Force as a column cell array of character vectors
yVars = force_column_cell(yVars);

% Extract the y vectors
yData = cellfun(@(x) dataTable{:, x}, yVars, 'UniformOutput', false);

% Creat full path to spreadsheet file if requested
if ~isempty(sheetName)
    sheetName = construct_fullpath(sheetName, 'Directory', outFolder);
end

%% Do the job
% Generate overlapped histograms
if plotHistograms
    % Create figure names if not provided
    if saveHistFlag && isempty(histFigNames)
        histFigNames = construct_fullpath(strcat('histogram_', yVars, '.png'), ...
                                        'Directory', outFolder);
    end

    % Plot grouped histograms
    h = cellfun(@(x, y, z) plot_grouped_histogram(x, xData, 'XLabel', y, ...
                        'FigName', z, 'Style', 'overlapped'), ...
                        yData, yVars, histFigNames, 'UniformOutput', false);
end

% Test for normality
% TODO

% Perform the appropriate comparison test
% TODO

%% Output results
% TODO: Place results in a table
% testResults = table(isDifferent, pValue);
testResults = table(histFigNames);
if plotHistograms
    testResults = addvars(testResults, histFigNames, h);
end

% Save the table if requested
if ~isempty(sheetName)
    writetable(testResults, sheetName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if isstring(yVars)
    % Use arrayfun
    yData = arrayfun(@(x) dataTable{:, x}, yVars, 'UniformOutput', false);
else
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%