function swpInfo = m3ha_load_sweep_info(varargin)
%% Loads sweep info (default is homeDirectory/data_dclamp/take4/dclampdatalog_take4.csv)
% Usage: swpInfo = m3ha_load_sweep_info(varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       swpInfo     - TODO: Description of swpInfo
%                   specified as a TODO
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
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
% TODO: Make homeDirectory an optional parameter

%% Hard-coded parameters
% Directories used
dataDirName = fullfile('data_dclamp', 'take4');
datalogFileName = 'dclampdatalog_take4.csv';

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, varargin{:});

%% Preparation
% Locate the home directory
homeDirectory = m3ha_locate_homedir;

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