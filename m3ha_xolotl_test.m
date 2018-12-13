% m3ha_xolotl_plot
%% Tests xolotl
%
% Requires:
%       cd/m3ha_xolotl_create_neuron.m
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_set_simparams.m
%       cd/m3ha_xolotl_plot.m

% File History:
% 2018-12-12 Created by Adam Lu

% Hard-coded parameters
neuronParamsFile = fullfile('/media/adamX/m3ha/optimizer4gabab', ...
                            'initial_params/initial_params_D091710.csv');
cpDelay = 1100;        % current pulse delay in ms
cpDuration = 10;       % current pulse duration in ms
cpAmplitude = -0.050;  % current pulse amplitude in nA

%{
timeStepDefault = 0.1;      % sample at 10 kHz by default
timeEndDefault = 1350;      % simulate for 1350 ms by default
simTimeStepDefault = [];    % same as timeStep by default
closedLoopDefault = false;  % don't use the final condition as the initial
                            %   condition for the next simulation by default
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a xolotl object based on a parameters file
m3ha = m3ha_xolotl_create_neuron(neuronParamsFile);

% Set simulation parameters
% TODO: set params here
m3ha = xolotl_set_simparams(m3ha);

% Add a current clamp
m3ha = xolotl_add_current_pulse(m3ha, 'Delay', cpDelay, ...
                            'Duration', cpDuration, 'Amplitude', cpAmplitude);

% Simulate and plot
%m3ha.plot;

% Simulate only
% vVecs = m3ha.integrate;

% Plot
% m3ha = m3ha_xolotl_plot(m3ha, vVecs);

% Save the xolotl object
% xolotl_save(m3ha, fileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%