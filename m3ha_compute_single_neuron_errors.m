function errors = m3ha_compute_single_neuron_errors (vSim, vReal, varargin)
%% Computes all errors for single neuron data
% Usage: errors = m3ha_compute_single_neuron_errors (vSim, vReal, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:    
%       vSim        - simulated voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       vReal       - recorded voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'TimeVecs': TODO: Description of TimeVecs
%                   must be a TODO
%                   default == TODO
%                   - 'IvecsSim': TODO: Description of IvecsSim
%                   must be a TODO
%                   default == TODO
%                   - 'IvecsReal': TODO: Description of IvecsReal
%                   must be a TODO
%                   default == TODO
%                   - 'SimMode': TODO: Description of SimMode
%                   must be a TODO
%                   default == TODO
%                   - 'FitWindow': TODO: Description of FitWindow
%                   must be a TODO
%                   default == TODO
%                   - 'BaseWindow': TODO: Description of BaseWindow
%                   must be a TODO
%                   default == TODO
%                   - 'BaseNoise': TODO: Description of BaseNoise
%                   must be a TODO
%                   default == TODO
%                   - 'SweepWeights': TODO: Description of SweepWeights
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/iscellnumeric.m
%       cd/match_row_counts.m
%
% Used by:    
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

% File History:
% 2018-10-24 Adapted from code in run_neuron_once_4compgabab.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vSim', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vSim must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'vReal', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vReal must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, vReal, varargin{:});
param1 = iP.Results.param1;

% Make sure windows are row vectors
if numel(fitWindow) == 2 && iscolumn(fitWindow)
    fitWindow = transpose(fitWindow);
end

%% Preparation
% Count the number of sweeps
nSweeps = numel(vSim);

% Simulation mode specific operations
switch simMode
    case 'passive'
    case 'active'
    otherwise
        error('simMode unrecognized!');
end

% Get the lower and upper bounds of the fit windows
fitWindow = match_row_counts(fitWindow, nSweeps);
fitWindowLowerBounds = fitWindow(:, 1);
fitWindowUpperBounds = fitWindow(:, 2);

% Calculate constants for efficiency
totalSweepWeights = sum(sweepWeights);
totalltsw = sum(ltsWeights);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fitWindowLowerBounds = match_vector_counts(fitWindowLowerBounds, [nSweeps, 1]);
fitWindowUpperBounds = match_vector_counts(fitWindowUpperBounds, [nSweeps, 1]);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%