% function m3ha_network_compare_ipsc (varargin)
%% Compare an evoked IPSC against the recorded IPSC
% Usage: m3ha_network_compare_ipsc (varargin)
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
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/compute_gabab_conductance.m
%       cd/create_labels_from_numbers.m
%       cd/extract_columns.m
%       cd/find_matching_files.m
%       cd/load_neuron_outputs.m
%       cd/m3ha_load_gabab_ipsc_params.m
%       cd/plot_traces.m
%       cd/set_figure_properties.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2020-01-22 Created by Adam Lu
% 

%% Hard-coded parameters
spExtension = 'singsp';
spPrefix = 'TC[0]';
spKeyword = 'gIncr_20';
ipscStartMs = 3000;
ampScaleFactor = 200;
ampUnits = 'nS';

% Column numbers for simulated data
%   Note: Must be consistent with m3ha_net.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
INA_COL_SIM = 3;
IK_COL_SIM = 4;
ICA_COL_SIM = 5;
IGABAA_COL_SIM = 6;
IGABAB_COL_SIM = 7;
CAI_COL_SIM = 8;
GGABAB_COL_SIM = 9;

% Plot parameters
xLimits = [2000, 10000];
xLabel = 'Time (ms)';
pharmLabels = {'{\it s}-Control', '{\it s}-GAT1 Block', ...
                    '{\it s}-GAT3 Block', '{\it s}-Dual Block'};

% TODO: Make optional arguments
directory = pwd;

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
% iP = inputParser;
% iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
% parse(iP, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Preparation
% TODO

%% Do the job
% Construct the pharm strings expected in file names
pharmStrs = create_labels_from_numbers(1:4, 'Prefix', 'pCond_');

% Locate the TC neuron data for each pharm condition
[~, dataPaths] = find_matching_files(pharmStrs, 'Directory', directory, ...
                            'Prefix', spPrefix, 'Keyword', spKeyword, ...
                            'Extension', spExtension);

% Load simulated data
simData = load_neuron_outputs('FileNames', dataPaths);

% Extract vectors from simulated data
[tVecsMs, gCmdSimUs] = extract_columns(simData, [TIME_COL_SIM, GGABAB_COL_SIM]);

% Convert to nS
gCmdSimNs = convert_units(gCmdSimUs, 'uS', 'nS');

% Load default GABAB IPSC parameters in nS
[ampOrig, tauRiseOrig, tauFallFastOrig, tauFallSlowOrig, weightOrig] = ...
    m3ha_load_gabab_ipsc_params('AmpScaleFactor', ampScaleFactor, ...
                                'AmpUnits', ampUnits);

% Compute original GABAB conductance vectors in nS
gVecsOrigNs = compute_gabab_conductance(tVecsMs, ipscStartMs, ...
                                ampOrig, tauRiseOrig, ...
                                tauFallFastOrig, tauFallSlowOrig, weightOrig);

% Clear simData to release memory
clear simData

% Create a figure
fig = set_figure_properties('AlwaysNew', true);

% Plot traces
plot_traces(tVecsMs, gCmdSimNs, 'DataToCompare', gVecsOrigNs, ...
            'PlotMode', 'parallel', ...
            'XLimits', xLimits, 'XLabel', xLabel, 'TraceLabels', pharmLabels);

%% Output results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%