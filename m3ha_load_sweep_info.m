function [swpInfo, varargout] = m3ha_load_sweep_info (varargin)
%% Loads sweep info (default is m3ha_locate_homedir/data_dclamp/take4/dclampdatalog_take4.csv)
% Usage: [swpInfo, cellInfo] = m3ha_load_sweep_info (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       swpInfo = m3ha_load_sweep_info;
%       [swpInfo, cellInfo] = m3ha_load_sweep_info;
%
% Outputs:
%       swpInfo     - a table containing sweep information
%                   specified as a table
%       cellInfo    - a table of cell info
%                   specified as a 2D table with fields specified in
%                       m3ha_create_cell_info_table.m
%
% Arguments:
%       varargin    - 'Directory': home directory containing sweep info
%                   must be a string scalar or a character vector
%                   default == m3ha_locate_homedir/data_dclamp/take4;
%                   - 'FileName': file name for the sweep info spreadsheet
%                   must be a string scalar or a character vector
%                   default == dclampdatalog_take4.csv
%
% Requires:
%       cd/m3ha_locate_homedir.m
%
% Used by:
%       cd/m3ha_compute_and_plot_statistics.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_create_cell_info_table.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_organize_sweep_indices.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_select_cells.m
%       cd/m3ha_select_raw_traces.m
%       cd/m3ha_select_sweeps.m
%       cd/m3ha_simulate_population.m

% File History:
% 2018-12-05 Adapted from code in singleneuronfitting42.m
% 2019-11-14 Made 'Directory' an optional parameter
% 2019-11-26 Made 'FileName' an optional parameter
% 2019-11-26 Now uses file bases as row names

%% Hard-coded parameters
% Directories used
dataDirName = fullfile('data_dclamp', 'take4');
datalogFileName = 'dclampdatalog_take4.csv';

%% Default values for optional arguments
directoryDefault = '';          % set later
fileNameDefault = '';          % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileName', fileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
fileName = iP.Results.FileName;

%% Preparation
% Decide on the home directory
if isempty(directory)
    directory = fullfile(m3ha_locate_homedir, dataDirName);
end

if isempty(fileName)
    fileName = datalogFileName;
end

% Decide on the full path to the data file
[dataPath, pathExists] = ...
    construct_and_check_fullpath(fileName, 'Directory', directory, ...
                                    'Extension', '.csv');

% Return if path doesn't exist
if ~pathExists
    fprintf('Warning: %s doesn''t exist!', dataPath);
    swpInfo = table.empty;
    return
end

%% Do the job
% Load sweep information from datalogPath
%   Note: the file names are read in as row names
swpInfo = readtable(dataPath, 'HeaderLines', 1, 'ReadVariableNames', true);

% Extract the file names
fileNames = swpInfo.fnrow;

% Extract the file bases
fileBases = extractBefore(fileNames, '.');

% Set the file bases as row names
swpInfo.Properties.RowNames = fileBases;

%% Create cell info table if requested
if nargout >= 2
    % Create cell info table
    cellInfo = m3ha_create_cell_info_table('SwpInfo', swpInfo);

    % Output
    varargout{1} = cellInfo;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
