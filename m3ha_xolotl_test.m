% m3ha_xolotl_test
%% Tests xolotl
%
% Requires:
% TODO
%       cd/extract_columns.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_xolotl_create_neuron.m
%       cd/m3ha_xolotl_plot.m
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_add_holding_current.m
%       cd/xolotl_set_simparams.m

% File History:
% 2018-12-12 Created by Adam Lu
% 2018-12-13 Added voltage clamp
% 2019-01-01 Now uses m3ha_import_raw_traces.m

%% Hard-coded parameters
% Files
sweepName = 'D091710_0012_15';
cellName = sweepName(1:7);
projectDir = '/media/adamX/m3ha';
realDataDir = fullfile(projectDir, 'data_dclamp');
matFilesDir = fullfile(realDataDir, 'take4/matfiles');
simDir = fullfile(projectDir, 'optimizer4gabab');
outFolder = '/media/adamX/xolotl_test';
outFileName = 'xolotl_test_before_simulations.mat';
neuronParamsDir = fullfile(simDir, 'initial_params');
neuronParamsFileName = ['initial_params_', cellName, '.csv'];
figTitle = ['Simulation for ', replace(sweepName, '_', '\_')];

% Parameters that should be consistent with the experiment
temperature = 33;       % temperature of 33 degrees Celsius used by Christine
compToPatch = 'soma';   % compartment to be patched
cprWindowOrig = [0, 350];   % original current pulse window in ms
cpStartWindowOrig = [95, 105];
cpDelayOrig = 100;      % original current pulse delay in ms
cpDuration = 10;        % current pulse duration in ms
cpAmplitude = -0.050;   % current pulse amplitude in nA

% Parameters for simulations
nCompartments = 3;      % number of compartments
closedLoop = false;     % whether to use the final state as the initial
                        %   condition for the next simulation
solverOrder = 0;        % uses the implicit Crank-Nicholson scheme
                        %   for multi-compartment models
timeStep = 0.1;         % output time step in ms
simTimeStep = 0.1;      % simulation time step in ms
timeToStabilize = 1000; % time for state variables to stabilize in ms

% Parameters for calculations
baseWindowOrig = [0, cpDelayOrig];
fitWindowOrig = [cpDelayOrig, cprWindowOrig(2)];
                        % from current pulse start to current pulse response end

% Parameters for plotting
xLimitsOrig = [cpDelayOrig - 10, cpDelayOrig + 40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Generate full file paths
neuronParamsPath = fullfile(neuronParamsDir, neuronParamsFileName);
outPath = fullfile(outFolder, outFileName);

% Shift original times by timeToStabilize
cprWindow = cprWindowOrig + timeToStabilize;
cpDelay = cpDelayOrig + timeToStabilize;
fitWindow = fitWindowOrig + timeToStabilize;
baseWindow = baseWindowOrig + timeToStabilize;
xLimits = xLimitsOrig + timeToStabilize;

% Obtain the simulation end time in ms for current pulse response
timeEndCpr = cprWindow(2);

%% Import data to compare against
% TODO: Use m3ha_import_raw_traces.m
[dataCpr, sweepInfo] = m3ha_import_raw_traces(sweepName, 'Directory', matFilesDir, ...
                        'Verbose', true, 'CreateLog', true, ...
                        'ToParsePulse', true, 'ToCorrectDcSteps', true, ...
                        'ToAverageByVhold', true, 'OutFolder', outFolder, ...
                        'TimeToPad', timeToStabilize, ...
                        'ResponseWindow', cprWindowOrig, ...
                        'StimStartWindow', cpStartWindowOrig);

% Extract needed data and information
[realVoltageData, realStimPulse] = extract_columns(dataCpr, 2:3);
holdingPotential = sweepInfo.holdPotential;

% Put together as data to compare
dataToCompare = [realVoltageData, realVoltageData, ...
                 realVoltageData, realStimPulse];

% Make sure holdingPotential has the correct dimensions
holdingPotential = match_dimensions(holdingPotential, [1, nCompartments]);

%% Create the neuron
% Create a xolotl object based on a parameters file
%   Note: m3ha is a handle to the xolotl object
m3ha = m3ha_xolotl_create_neuron(neuronParamsPath);

% Set general simulation parameters
xolotl_set_simparams(m3ha, 'ClosedLoop', closedLoop, ...
                    'SolverOrder', solverOrder, 'Temperature', temperature, ...
                    'TimeStep', timeStep, 'SimTimeStep', simTimeStep);

% Set simulation parameters for the current pulse response
xolotl_set_simparams(m3ha, 'TimeEnd', timeEndCpr, ...
                        'InitialVoltage', holdingPotential);

% Add a current pulse protocol
%   Note: must do this after the correct time_end is set
xolotl_add_current_pulse(m3ha, 'Compartment', compToPatch, ...
                                'Delay', cpDelay, 'Duration', cpDuration, ...
                                'Amplitude', cpAmplitude);

% Save the xolotl object before simulations
save(outPath, 'm3ha');

%% Simulate and plot
% Simulate and use default plot
% m3ha.plot;

% Simulate and plot individual traces against data
m3ha_xolotl_plot(m3ha, 'DataToCompare', dataToCompare, ...
                        'XLimits', xLimits, ...
                        'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
                        'FigTitle', figTitle, 'CompToPatch', compToPatch, ...
                        'TimeToStabilize', timeToStabilize, ...
                        'HoldingPotential', holdingPotential);


% Set up the plots to be manipulated
m3ha.manipulate_plot_func = {@m3ha_xolotl_plot};

% Manipulate leak channel parameters
m3ha.manipulate('*gbar')
% m3ha.manipulate('*Leak*')
% m3ha.manipulate('*E')
% m3ha.manipulate({'*Leak*', '*length*'})
manip = gcf;
manip.Position(2) = manip.Position(2) - 200;

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

% TODO: xolotl_save(xolotlObject, 'OutFolder', outFolder, 'FileBase', fileBase);

holdingPotential = -60; % holding potential in mV

% Compute baseline window
if isempty(baseWindowOrig)
    baseWindowOrig = compute_default_sweep_info(tvecOrig, vvecCpr, ...
                                            'FitWindow', fitWindowOrig);
end

matFile = fullfile(matFilesDir, [sweepName, '.mat']);

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
endPointsCpr = find_window_endpoints(cprWindowOrig, tvecOrig);

% Extract the current pulse response window
[tvecCpr, vvecCpr, ivecCpr] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsCpr), ...
            tvecOrig, vvecOrig, ivecOrig);

% Find the end points for the baseline window
endPointsBase = find_window_endpoints(baseWindowOrig, tvecOrig);
                                    
% Compute the holding potential
holdingPotential = compute_means(vvecOrig, 'EndPoints', endPointsBase);

% Repeat the current pulse responses for each compartment for now
vvecsCpr = repmat(vvecCpr, [1, nCompartments]);

% Construct padding ones or zeros
vvecsToPadCpr = ones(nSamplesToPadForStabilization, 1) * holdingPotential;
ivecToPadCpr = zeros(nSamplesToPadForStabilization, 1);

realVoltageData = [vvecsToPadCpr; vvecsCpr];
realStimPulse = [ivecToPadCpr; ivecCpr];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%