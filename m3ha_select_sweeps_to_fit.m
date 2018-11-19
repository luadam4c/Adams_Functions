function [fileNamesToFit, swpIndicesToFit] = ...
                m3ha_select_sweeps_to_fit (swpInfo, varargin)
%% Find file names and row indices in swpInfo that will be used for fitting
% Usage: [fileNamesToFit, swpIndicesToFit] = ...
%               m3ha_select_sweeps_to_fit (swpInfo, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       fileNamesToFit  - file names of sweeps to fit
%                       specified as a cell array
%       swpIndicesToFit - row indices in swpInfo of sweeps to fit
%                       specified as a positive integer vector
% Arguments:
%       swpInfo     - a table of sweep info, with each row named by 
%                       the matfile name containing the raw data
%                   must be a 2D table with row names being file names
%                       and with the fields:
%                       cellidrow - cell ID
%                       prow      - pharmacological condition
%                       grow      - conductance amplitude scaling (%)
%       varargin    - 'DataMode': data mode
%                   must be a one of:
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                           but exclude cell-pharm-g_incr sets 
%                           containing problematic sweeps
%                   default == 2
%                   - 'CasesDir' - the directory that contains 
%                                   'TAKE_OUT_*' folders with special cases
%                   must be a directory
%                   default == ~/m3ha/data_dclamp/take4/special_cases
%
% Requires:
%       cd/find_rows_with_same_attributes.m
%       cd/m3ha_find_files_to_take_out.m
%
% Used by:
%       ~/m3ha/optimizer4gabab/singleneuronfitting42.m 

% TODO:
%       /media/adamX/m3ha/data_dclamp/dclampPassiveFitter.m
%       /media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%       /media/adamX/m3ha/data_dclamp/PlotCorrelations.m
%       /media/adamX/m3ha/data_dclamp/dclampdatalog_analyze.m

% File History:
% 2018-11-18 Adapted from m3ha_find_ind_to_fit()
% 

%% Hard-coded parameters
attributesToMatch = {'cellidrow', 'prow', 'grow'};

%% Default values for optional arguments
dataModeDefault = 2;
casesDirDefault = '~/m3ha/data_dclamp/take4/special_cases';

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
addRequired(iP, 'swpInfo', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));
addParameter(iP, 'CasesDir', casesDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, swpInfo, varargin{:});
dataMode = iP.Results.DataMode;
casesDir = iP.Results.CasesDir;

%% Preparation
% Extract all the row names
fileNames = swpInfo.Properties.RowNames;

% Extract the conductance amplitude scaling percentages
gIncr = swpInfo.grow;

%% Do the job
% Find the sweep indices for 
%   conductance amplitudes with 100%, 200% or 400% scaling 
%   (these are present in all experiments)
swpIndGIncrToFit = find(gIncr == 100 | gIncr == 200 | gIncr == 400);

% Get the file names of files to take out from specialCasesDir
%   Note: these are labelled with 'TAKE_OUT_*' and were
%           the result of voting by blab
if dataMode == 2
    % Find all files to take out
    fileNamesToTakeOut = m3ha_find_files_to_take_out('CasesDir', casesDir);

    % If there are no files, return with message
    if isempty(fileNamesToTakeOut)
        fprintf('There are no files to take out under %s!\n', casesDir);
        fileNamesToFit = {};
        swpIndicesToFit = [];
        return
    end
        
    % Find all the sweep indices with the same cell ID, pharm condition
    %   and conductance amplitude scaling as the files to take out
    swpIndNotToFit = ...
        find_rows_with_same_attributes(swpInfo, fileNamesToTakeOut, ...
                                        attributesToMatch);
end

% Find the sweep indices to fit
if dataMode == 1
    swpIndicesToFit = swpIndGIncrToFit;
elseif dataMode == 2
    swpIndicesToFit = setdiff(swpIndGIncrToFit, swpIndNotToFit);
end

% Find the file names to fit
fileNamesToFit = fileNames(swpIndicesToFit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find the indices of all sweeps of the same cell-pharm-g incr group
indSameCellId = find(cellId == cellId(iSwp));
indSamePharm = find(pharm == pharm(iSwp));
indSameGIncr = find(gIncr == gIncr(iSwp));
indSameGroup = intersect(indSameCellId, ...
                intersect(indSamePharm, indSameGIncr));

swpIndSameGroup = cell(nSweeps, 1);

parfor iSwp = 1:nSweeps
    if ismember(fileNames(iSwp), fileNamesToTakeOut)
        % If a sweep to take out, find the indices of all sweeps 
        %   of the same cell-pharm-g incr group
        swpIndSameGroup{iSwp} = ...
            find(cellId == cellId(iSwp) & pharm == pharm(iSwp) & ...
                gIncr == gIncr(iSwp));
    else
        % Otherwise, return nothing
        swpIndSameGroup{iSwp} = [];
    end
end

swpIndNotToFit = unique(union_over_cells(swpIndSameGroup));

% Check if this sweep is not yet included in swpIndNotToFit
isNotIncluded = isempty(find(swpIndNotToFit == iSwp, 1));

% Sort swpIndNotToFit
swpIndNotToFit = sort(swpIndNotToFit);

% Extract from table
cellId = table.cellidrow;
pharm = table.prow;
gIncr = table.grow;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%