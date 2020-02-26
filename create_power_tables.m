function [sampsizeByPowerTable, sampsizeByEffectTable] = ...
                create_power_tables (keyword, directory, statPowerDesired, ...
                                    useLog2Ratio, meanAltExpected, ...
                                    measureString, varName1, varName2, varargin)
%% Creates power tables
% Usage: [sampsizeByPowerTable, sampsizeByEffectTable] = ...
%               create_power_tables (keyword, directory, statPowerDesired, ...
%                                   useLog2Ratio, meanAltExpected, ...
%                                   measureString, varName1, varName2, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for compute_sampsizepwr()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/compute_sampsizepwr.m
%       cd/find_matching_files.m
%
% Used by:
%       cd/metabolismR01_power_analysis.m

% File History:
% 2020-02-25 Moved from plethR01_power_analysis.m
% 

%% Hard-coded parameters
% Output file names
byPowerSuffix = 'power_table_by_power';
byEffectSuffix = 'power_table_by_effect_size';

% TODO: Make optional arguments
sheetType = 'csv';

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 8
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'keyword');
addRequired(iP, 'directory');
addRequired(iP, 'statPowerDesired');
addRequired(iP, 'useLog2Ratio');
addRequired(iP, 'meanAltExpected');
addRequired(iP, 'measureString');
addRequired(iP, 'varName1');
addRequired(iP, 'varName2');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, keyword, directory, statPowerDesired, ...
            useLog2Ratio, meanAltExpected, ...
            measureString, varName1, varName2, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the compute_sampsizepwr() function
otherArguments = iP.Unmatched;

%% Preparation
% Find the spreadsheet of interest
% TODO: Improve find_matching_files to take 'Extension' as a cell array
[~, dataPath] = find_matching_files(keyword, 'Directory', directory, ...
                                    'Extension', 'csv', 'PartType', 'suffix');
if isempty(dataPath)
    [~, dataPath] = find_matching_files(keyword, 'Directory', directory, ...
                                    'Extension', 'xlsx', 'PartType', 'suffix');
end
if isempty(dataPath)
    [~, dataPath] = find_matching_files(keyword, 'Directory', directory, ...
                                    'Extension', 'xls', 'PartType', 'suffix');
end

% Set log ratio info
if useLog2Ratio
    measureString = [measureString, '_log2ratio'];
end

% Create paths
[byPowerPath, byEffectPath] = ...
    argfun(@(x) fullfile(directory, [keyword, '_', measureString, ...
                                    '_', x, '.', sheetType]), ...
            byPowerSuffix, byEffectSuffix);

% Read in data
dataTable = readtable(dataPath);
data1 = dataTable.(varName1);
data2 = dataTable.(varName2);

% Compute log2 data
if useLog2Ratio
    log2Data1Orig = log2(data1);
    log2Data2Orig = log2(data2);

    % Remove any infinite data
    % TODO: Move this to compute_sampsizepwr.m?
    data1 = log2Data1Orig(~isinf(log2Data1Orig) & ~isinf(log2Data2Orig));
    data2 = log2Data2Orig(~isinf(log2Data1Orig) & ~isinf(log2Data2Orig));
end

%% Do the job
% Compute the power table by varying statistical power
sampsizeByPowerTable = ...
    compute_sampsizepwr('DataCon', data1, 'DataExp', data2, ...
            'StatPowerDesired', statPowerDesired, 'OutputType', 'scalars', ...
            otherArguments);

% Compute the power table by varying effect size
sampsizeByEffectTable = ...
    compute_sampsizepwr('DataCon', data1, 'DataExp', data2, ...
            'MeanAlt', meanAltExpected, 'OutputType', 'scalars', ...
            otherArguments);

%% Save the tables
writetable(sampsizeByPowerTable, byPowerPath);
writetable(sampsizeByEffectTable, byEffectPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%