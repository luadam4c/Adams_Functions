function handles = m3ha_compute_and_plot_violin (directory, ...
                                pharmLabels, conditionLabel, varargin)
%% Computes the statistics for and plots 2D violin plots
% Usage: handles = m3ha_compute_and_plot_violin (statsPath, directory, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - structures with fields:
%                       fig
%                       violins
%                   specified as a structure array
%
% Arguments:
%       directory       - TODO: Description of statsPath
%                       must be a TODO
%       pharmLabels     - TODO: Description of statsPath
%                       must be a TODO
%       conditionLabel  - TODO: Description of statsPath
%                       must be a TODO
%       varargin    - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow       - cell ID
%                       ltspeaktimes    - LTS peak delays
%                       spikesperpeak   - action potentials per LTS
%                       bursttimes      - burst delays
%                       spikesperburst  - action potentials per burst
%                   default == m3ha_load_sweep_info
%                   - 'DataMode': data mode
%                   must be one of:
%                       0 - all data
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == 0
%                   - 'PharmConditions': pharmacological condition(s)
%                                           to restrict to
%                   must be empty or some of:
%                       1 - control
%                       2 - GAT1 blockade
%                       3 - GAT3 blockade
%                       4 - dual blockade
%                       or a cell array of them (will become 1st dimension)
%                   default == no restrictions
%                   - 'GIncrCondition': conductance amplitude condition(s) (%)
%                                           to restrict to
%                   must be empty or some of: 25, 50, 100, 200, 400, 800
%                       or a cell array of them (will become 2nd dimension)
%                   default == no restrictions
%                   - 'VHoldConditions': holding potential condition(s) (mV)
%                                           to restrict to
%                   must be empty or some of: -60, -65, -70
%                       or a cell array of them (will become 3rd dimension)
%                   default == no restrictions
%                   - Any other parameter-value pair for m3ha_plot_violin()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_plot_violin.m
%
% Used by:
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_simulate_population.m

% File History:
% 2020-03-11 Adapted from m3ha_plot_figure02.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
swpInfoDefault = table.empty;
dataModeDefault = 0;
pharmConditionsDefault = [];
gIncrConditionsDefault = [];
vHoldConditionsDefault = [];

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
addRequired(iP, 'directory');
addRequired(iP, 'pharmLabels');
addRequired(iP, 'conditionLabel');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
addParameter(iP, 'PharmConditions', pharmConditionsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'GIncrConditions', gIncrConditionsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'VHoldConditions', vHoldConditionsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));

% Read from the Input Parser
parse(iP, directory, pharmLabels, conditionLabel, varargin{:});
swpInfo = iP.Results.SwpInfo;
dataMode = iP.Results.DataMode;
pharmConditions = iP.Results.PharmConditions;
gIncrConditions = iP.Results.GIncrConditions;
vHoldConditions = iP.Results.VHoldConditions;

% Keep unmatched arguments for the m3ha_plot_violin() function
otherArguments = iP.Unmatched;

%% Do the job
% Construct stats table path
statsPath = fullfile(directory, strcat(conditionLabel, '_stats.mat'));

% Compute statistics if not done already
if ~isfile(statsPath)
    % Compute statistics for all features
    disp('Computing statistics for violin plots ...');
    statsTable = m3ha_compute_statistics('SwpInfo', swpInfo, ...
                                        'DataMode', dataMode, ...
                                        'PharmConditions', pharmConditions, ...
                                        'GIncrConditions', gIncrConditions, ...
                                        'VHoldConditions', vHoldConditions);

    % Save stats table
    save(statsPath, 'statsTable', 'pharmLabels', ...
                        'conditionLabel', '-v7.3');
end

% Plot violin plots
m3ha_plot_violin(statsPath, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%