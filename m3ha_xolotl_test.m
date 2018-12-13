% m3ha_xolotl_test
%% Tests xolotl
%
% Requires:
%       cd/m3ha_xolotl_create_neuron.m
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_add_holding_current.m
%       cd/xolotl_add_voltage_clamp.m
%       cd/xolotl_set_simparams.m
%       cd/xolotl_simulate_voltage_clamp.m
%       cd/m3ha_xolotl_plot.m

% File History:
% 2018-12-12 Created by Adam Lu
% 2018-12-13 Added voltage clamp

% Hard-coded parameters
neuronParamsFile = fullfile('/media/adamX/m3ha/optimizer4gabab', ...
                            'initial_params/initial_params_D091710.csv');
cpDelay = 1100;         % current pulse delay in ms
cpDuration = 10;        % current pulse duration in ms
cpAmplitude = -0.050;   % current pulse amplitude in nA
timeStep = 0.1;         % output time step in ms
initialVoltageVc = -70; % initial voltage level for voltage clamp
initialVoltageCpr = -60;% initial voltage level for current pulse response
timeEndVc = 1000;       % simulation end time in ms for voltage clamp
timeEndCpr = 1350;      % simulation end time in ms for current pulse response
simTimeStep = 0.1;      % simulation time step in ms
closedLoop = false;     % whether to use the final state as the initial
                        %   condition for the next simulation
solverOrder = 0;        % uses the implicit Crank-Nicholson scheme
                        %   for multi-compartment models
temperature = 33;       % temperature of 33 degrees Celsius used by Christine
holdingPotential = -60; % holding potential in mV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Create a xolotl object based on a parameters file
m3ha = m3ha_xolotl_create_neuron(neuronParamsFile);

% Set simulation parameters
m3ha = xolotl_set_simparams(m3ha, 'TimeStep', timeStep, ...
                    'SimTimeStep', simTimeStep, 'ClosedLoop', closedLoop, ...
                    'Temperature', temperature, 'SolverOrder', solverOrder);

%% Find the holding current
% Set simulation parameters for voltage clamp
m3ha = xolotl_set_simparams(m3ha, 'InitialVoltage', initialVoltageVc, ...
                            'TimeEnd', timeEndVc);

% Add a voltage clamp
m3haTest = xolotl_add_voltage_clamp(m3ha, 'Amplitude', holdingPotential);

% Find the holding current necessary to match the holding potential
% TODO: holdingCurrent = xolotl_simulate_voltage_clamp(m3haTest);
holdingCurrent = 0;

%% Simulate the current pulse protocol
% Set simulation parameters for current pulse response
m3ha = xolotl_set_simparams(m3ha, 'InitialVoltage', initialVoltageCpr, ...
                            'TimeEnd', timeEndCpr);

% Add a holding current
m3ha = xolotl_add_holding_current(m3ha, 'Amplitude', holdingCurrent);

% Add a current pulse protocol
m3ha = xolotl_add_current_pulse(m3ha, 'Delay', cpDelay, ...
                            'Duration', cpDuration, 'Amplitude', cpAmplitude);

% Simulate and plot
% m3ha.plot;

% Simulate only
% vVecs = m3ha.integrate;

% Plot
% m3ha = m3ha_xolotl_plot(m3ha, vVecs);

% Save the xolotl object
% xolotl_save(m3ha, fileName);

% Displays a list of properties
% properties(xolotl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%