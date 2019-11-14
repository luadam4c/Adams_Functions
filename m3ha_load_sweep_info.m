function swpInfo = m3ha_load_sweep_info(varargin)
%% Loads sweep info (default is homeDirectory/data_dclamp/take4/dclampdatalog_take4.csv)
% Usage: swpInfo = m3ha_load_sweep_info(varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       swpInfo     - a table containing sweep information
%                   specified as a table
% Arguments:
%       varargin    - 'HomeDirectory': home directory containing sweep info
%                   must be a string scalar or a character vector
%                   default == m3ha_locate_homedir;
%
% Requires:
%       cd/m3ha_locate_homedir.m
%
% Used by:
%       cd/m3ha_generate_cell_info.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_organize_sweep_indices.m
%       cd/m3ha_select_cells.m
%       cd/m3ha_select_sweeps_to_fit.m

% File History:
% 2018-12-05 Adapted from code in singleneuronfitting42.m
% 2019-11-14 Made homeDirectory an optional parameter

%% Hard-coded parameters
% Directories used
dataDirName = fullfile('data_dclamp', 'take4');
datalogFileName = 'dclampdatalog_take4.csv';

%% Default values for optional arguments
homeDirectoryDefault = [];      % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'HomeDirectory', homeDirectoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
homeDirectory = iP.Results.HomeDirectory;

%% Preparation
% Decide on the home directory
if isempty(homeDirectory)
    homeDirectory = m3ha_locate_homedir;
end

% Generate full path to data
dataDir = fullfile(homeDirectory, dataDirName);
datalogPath = fullfile(dataDir, datalogFileName);

%% Do the job
% Load sweep information from datalogPath
%   Note: the file names are read in as row names
swpInfo = readtable(datalogPath, 'HeaderLines', 1, ...
                    'ReadVariableNames', true, 'ReadRowNames', true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%