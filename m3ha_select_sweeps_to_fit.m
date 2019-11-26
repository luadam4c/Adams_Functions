function [swpInfo, fileBasesToFit] = m3ha_select_sweeps_to_fit (varargin)
%% Find file names and row indices in swpInfo that will be used for fitting
% Usage: [swpInfo, fileBasesToFit] = m3ha_select_sweeps_to_fit (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       swpInfo         - same as input swpInfo table but with the additional field:
%                           toFit   - whether the sweep is to be fitted
%       fileBasesToFit  - file bases of sweeps to fit
%                       specified as a cell array
% Arguments:
%       varargin    - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                   default == m3ha_load_sweep_info
%                   - 'DataMode': data mode
%                   must be a one of:
%                       0 - all data
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == 2
%                   - 'CasesDir' - the directory that contains 
%                                   'TAKE_OUT_*' folders with special cases
%                   must be a directory
%                   default == ~/m3ha/data_dclamp/take4/special_cases
%
% Requires:
%       cd/has_same_attributes.m
%       cd/m3ha_find_files_to_take_out.m
%       cd/m3ha_load_sweep_info.m
%
% Used by:
%       cd/m3ha_compute_and_plot_statistics.m
%       cd/m3ha_compute_ltsburst_statistics.m
%       cd/m3ha_select_cells.m

% TODO:
%       /media/adamX/m3ha/data_dclamp/dclampPassiveFitter.m
%       /media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%       /media/adamX/m3ha/data_dclamp/PlotCorrelations.m
%       /media/adamX/m3ha/data_dclamp/dclampdatalog_analyze.m

% File History:
% 2018-11-18 Adapted from m3ha_find_ind_to_fit()
% 2018-12-06 Now adds a toFit column to sweep info
% 2019-11-26 Added dataMode == 0

%% Hard-coded parameters
attributesToMatch = {'cellidrow', 'prow', 'grow'};

%% Default values for optional arguments
swpInfoDefault = [];
dataModeDefault = 2;
casesDirDefault = '~/m3ha/data_dclamp/take4/special_cases';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));
addParameter(iP, 'CasesDir', casesDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, varargin{:});
swpInfo = iP.Results.SwpInfo;
dataMode = iP.Results.DataMode;
casesDir = iP.Results.CasesDir;

%% Preparation
% Read in swpInfo if not provided
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

% Extract the conductance amplitude scaling percentages
grow = swpInfo.grow;

%% Do the job
% Print message
fprintf('Selecting the sweeps to fit ... \n');

% Determine whether each sweep has
%   conductance amplitudes with 100%, 200% or 400% scaling 
%   (these are present in all experiments)
isGcondToFit = grow == 100 | grow == 200 | grow == 400;

% Get the file names of files to take out from specialCasesDir
%   Note: these are labelled with 'TAKE_OUT_*' and were
%           the result of voting by blab
if dataMode == 2
    % Find all files to take out
    fileNamesToTakeOut = m3ha_find_files_to_take_out('CasesDir', casesDir);

    % If there are no files, return with message
    if isempty(fileNamesToTakeOut)
        fprintf('There are no files to take out under %s!\n', casesDir);
        fileBasesToFit = {};
        swpIndicesToFit = [];
        return
    end

    % Determine whether each sweep as the same cell ID, pharm condition
    %   and conductance amplitude scaling as any of the files to take out
    isNotToFit = has_same_attributes(swpInfo, fileNamesToTakeOut, ...
                                        attributesToMatch);
end

% Find the sweep indices to fit
if dataMode == 0
    toFit = true(height(swpInfo), 1);
elseif dataMode == 1
    toFit = isGcondToFit;
elseif dataMode == 2
    toFit = isGcondToFit & ~isNotToFit;
end

%% Output results
% Place in swpInfo
swpInfo = addvars(swpInfo, toFit);

% Extract all the file bases
fileBases = swpInfo.Properties.RowNames;

% Find the file bases to fit
fileBasesToFit = fileBases(toFit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%