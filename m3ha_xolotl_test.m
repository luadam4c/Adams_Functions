% m3ha_xolotl_test
%% Tests xolotl
%
% Requires:
% TODO
%       cd/m3ha_xolotl_create_neuron.m
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_add_holding_current.m
%       cd/xolotl_set_simparams.m
%       cd/m3ha_xolotl_plot.m

% File History:
% 2018-12-12 Created by Adam Lu
% 2018-12-13 Added voltage clamp

% Hard-coded parameters
projectDir = '/media/adamX/m3ha';
realDataDir = fullfile(projectDir, 'data_dclamp');
matFilesDir = fullfile(realDataDir, 'take4/matfiles');
matFile = fullfile(matFilesDir, 'D091710_0012_15');
simDir = fullfile(projectDir, 'optimizer4gabab');
outFolder = '/media/adamX/xolotl_test';
outFileName = 'xolotl_test.mat';
neuronParamsDir = fullfile(simDir, 'initial_params');
neuronParamsFileName = 'initial_params_D091710.csv';
closedLoop = false;     % whether to use the final state as the initial
                        %   condition for the next simulation
solverOrder = 0;        % uses the implicit Crank-Nicholson scheme
                        %   for multi-compartment models
temperature = 33;       % temperature of 33 degrees Celsius used by Christine
timeStep = 0.1;         % output time step in ms
simTimeStep = 0.1;      % simulation time step in ms

compToPatch = 'soma';   % compartment to be patched
holdingPotential = -60; % holding potential in mV
timeToStabilize = 1000; % time for state variables to stabilize in ms

% TODO
cprWindow = [0, 350];
timeEndCpr = 1350;      % simulation end time in ms for current pulse response

cpDelay = 1100;         % current pulse delay in ms
cpDuration = 10;        % current pulse duration in ms
cpAmplitude = -0.050;   % current pulse amplitude in nA
xLimits = [1000, 1250];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Generate full file paths
neuronParamsPath = fullfile(neuronParamsDir, neuronParamsFileName);
outPath = fullfile(outFolder, outFileName);

%% Import data to compare against
% TODO: Make this its own function
%   [realVoltageData, realStimPulse] = m3ha_load_single_sweep(matFile)

% Open the matfile
m = matfile(matFile);

% Use the original, non-filtered traces
dataOrig = m.d_orig;

% Extract original data vectors
[tvecOrig, ivecOrig, vvecOrig] = extract_columns(dataOrig, [1, 3, 4]);

% Convert to nA
PA_PER_NA = 1000;
ivecOrig = ivecOrig/PA_PER_NA;

% Compute the sampling intervals in ms
siMsOrig = compute_sampling_interval(tvecOrig);

% Convert times to samples
[nSamplesToPadForStabilization] = ...
    argfun(@(x) convert_to_samples(x, siMsOrig), timeToStabilize);

% Find the time to stabilize
endPointsCpr = find_window_endpoints(cprWindow, tvecOrig);

% Extract the current pulse response window
[vvecCpr, ivecCpr] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsCpr), ...
            vvecOrig, ivecOrig);

% Construct 
vvecToPadCpr = ones(nSamplesToPadForStabilization, 1) * holdingPotential;
ivecToPadCpr = zeros(nSamplesToPadForStabilization, 1);

% Extract needed data
realVoltageData = [vvecToPadCpr; vvecCpr];
realStimPulse = [ivecToPadCpr; ivecCpr];

% Put together as data to compare
dataToCompare = [realVoltageData, realVoltageData, realVoltageData, realStimPulse];

%% Create the neuron
% Create a xolotl object based on a parameters file
m3ha = m3ha_xolotl_create_neuron(neuronParamsPath);

% Set general simulation parameters
m3ha = xolotl_set_simparams(m3ha, 'ClosedLoop', closedLoop, ...
                    'SolverOrder', solverOrder, 'Temperature', temperature, ...
                    'TimeStep', timeStep, 'SimTimeStep', simTimeStep);

%% Find the holding current (nA) necessary to match the holding potential
% holdingCurrent = ...
%     xolotl_estimate_holding_current(m3ha, holdingPotential, ...
%                                     'CompToPatch', compToPatch, ...
%                                     'TimeToStabilize', timeToStabilize);
% holdingCurrent = 0.0743;
holdingCurrent = 0.1125;
% holdingCurrent = 0;

%% Simulate the current pulse protocol
% Set simulation parameters for current pulse response
m3ha = xolotl_set_simparams(m3ha, 'InitialVoltage', holdingPotential, ...
                            'TimeEnd', timeEndCpr);

% Add a holding current (nA)
m3ha = xolotl_add_holding_current(m3ha, 'Compartment', compToPatch, ...
                                    'Amplitude', holdingCurrent);

% Add a current pulse protocol
m3ha = xolotl_add_current_pulse(m3ha, 'Compartment', compToPatch, ...
                                'Delay', cpDelay, 'Duration', cpDuration, ...
                                'Amplitude', cpAmplitude);

% Simulate and plot
% m3ha.plot;

% Simulate and plot against data
m3ha = m3ha_xolotl_plot(m3ha, 'DataToCompare', dataToCompare, ...
                        'XLimits', xLimits);


% Save the xolotl object
% TODO: xolotl_save(xolotlObject, 'OutFolder', outFolder, 'FileBase', fileBase);
save(outPath, 'm3ha');

mh3a.manipulate_plot_func = {@m3ha_xolotl_plot};

% Manipulate leak channel parameters
% m3ha.manipulate('*Leak*')
m3ha.manipulate('*gbar')
% m3ha.manipulate('*E')
% m3ha.manipulate({'*Leak*', '*length*'})

% Displays a list of properties
% properties(xolotl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

initialVoltageVc = -70; % initial voltage level for voltage clamp
timeEndVc = 1000;       % simulation end time in ms for voltage clamp

initialVoltageCpr = -60;% initial voltage level for current pulse response
m3ha = xolotl_set_simparams(m3ha, 'InitialVoltage', initialVoltageCpr, ...
                            'TimeEnd', timeEndCpr);

[holdingCurrent, testObject] = ...

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%