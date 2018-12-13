function xolotlObject = m3ha_xolotl_plot (xolotlObject, varargin)
%% Plots the simulation results from a xolotl object
% Usage: xolotlObject = m3ha_xolotl_plot (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       must be a xolotl object
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/count_samples.m
%       cd/create_time_vectors.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 201X-XX-XX Created by TODO or Adapted from TODO
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'xolotlObject', ...                  % TODO: Description of xolotlObject
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% TODO

% Get voltage traces for all compartments
vVecs = xolotlObject.integrate;

%% Plot results
% Get the number of samples
nSamples = count_samples(vVecs);

% Create a time vector in seconds
tVec = create_time_vectors(nSamples, 'SamplingIntervalMs', timeStep, ...
                            'TimeUnits', 's');

% Create figure
figure('Outerposition', [300, 300, 1200, 600], ...
        'PaperUnits', 'points', 'PaperSize', [1200, 600]);
hold on

% Create first subplot
subplot(3, 1, 1); hold on
plot(tVec, vVecs(:, 1), 'k')
xlabel('Time (s)')
ylabel('Voltage in soma (mV)')
set(gca, 'YLim', [-80 50])

% Create second subplot
subplot(3, 1, 2); hold on
plot(tVec, vVecs(:, 2),'k')
xlabel('Time (s)')
ylabel('Voltage in dendrite 1 (mV)')
set(gca, 'YLim', [-80 50])

% Create third subplot
subplot(3, 1, 3); hold on
plot(tVec, vVecs(:, 3),'k')
xlabel('Time (s)')
ylabel('Voltage in dendrite 1 (mV)')
set(gca, 'YLim', [-80 50])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

tVec = (1:length(vVecs)) * xolotlObject.dt * 1e-3;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%