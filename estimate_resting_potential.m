function [epas, Rin] = estimate_resting_potential (holdPotential, holdCurrentPa)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: [epas, Rin] = estimate_resting_potential (holdPotential, holdCurrentPa)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       /TODO:dir/TODO:file
%
% Used by:
%       cd/find_passive_params.m

% File History:
% 2018-11-13 Created by Adam Lu
% 

%% Hard-coded constants
PA_PER_NA = 1000;

%% Default values for optional arguments
param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'reqarg1', ...                  % TODO: Description of reqarg1
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, reqarg1, varargin{:});
param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% Force all vectors to be column vectors
holdCurrentNa = holdCurrentNa(:);
holdCurrentPa = holdCurrentPa(:);

%% Do the job
% Convert the holding current from pA to nA
holdCurrentNa = holdCurrentPa / PA_PER_NA;

% Estimate epas and Rin with linear least squares
% Define the matrices
% Ohm's Law: I * R + epas = V
%     units: [nA] * [MOhm] + [mV] = [mV]
% X * w = V
X = ones(numswps, 2);
X(:, 1) = holdCurrentNa;
V = holdPotential;

% Compute the estimates
w = pinv(X) * V;

%% Output results
% Extract the input resistance (MOhm)
Rin = w(1);         

% Extract the resting membrane potential (mV)
epas = w(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute the regression coefficient
RinRegression = (holdCurrent - mean(holdCurrent)) \ ...
                    (holdPotential - mean(holdPotential));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%