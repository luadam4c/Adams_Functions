% m3ha_xolotl_test
%% Tests xolotl
%
% Requires:
% TODO
%       cd/archive_dependent_scripts.m
%       cd/create_time_stamp.m
%       cd/compute_sampling_interval.m
%       cd/extract_columns.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_xolotl_create_neuron.m
%       cd/m3ha_xolotl_plot.m
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_add_holding_current.m
%       cd/xolotl_set_simparams.m
%       cd/xolotl_create_model_soplata.m
%       cd/xolotl_create_model_howard.m

% File History:
% 2018-12-12 Created by Adam Lu
% 2018-12-13 Added voltage clamp
% 2019-01-01 Now uses m3ha_import_raw_traces.m
% 2019-08-15 Now uses a time step that matches the sampling interval
% 2019-08-15 Now sets dendritic dataToCompare to be NaN values
% 2019-08-15 Added active mode
% 2019-08-15 Now uses current vector from data
% TODO: single compartment, try model_soplata

%% Hard-coded parameters
% Files
sweepName = 'C101210_0006_3';
% sweepName = 'D091710_0012_15';
cellName = sweepName(1:7);
projectDir = '/media/adamX/m3ha';
realDataDir = fullfile(projectDir, 'data_dclamp');
matFilesDir = fullfile(realDataDir, 'take4/matfiles');
simDir = fullfile(projectDir, 'optimizer4gabab');
outFolder = '/media/adamX/xolotl_test';
neuronParamsDir = fullfile(simDir, 'best_params', ...
            'bestparams_20180424_singleneuronfitting21_Rivanna');
% neuronParamsDir = fullfile(simDir, 'initial_params');
neuronParamsFileName = ['bestparams_', cellName, '.csv'];
% neuronParamsFileName = ['initial_params_', cellName, '.csv'];
figTitle = ['Simulation for ', replace(sweepName, '_', '\_')];

% Parameters that should be consistent with the experiment
temperature = 33;       % temperature of 33 degrees Celsius used by Christine
compToPatch = 'soma';   % compartment to be patched
cprWindowOrig = [0, 350];   % original current pulse window in ms
cpStartWindowOrig = [95, 105];
cpDelayOrig = 100;      % original current pulse delay in ms
cpDuration = 10;        % current pulse duration in ms
cpAmplitude = -0.050;   % current pulse amplitude in nA
ipscrWindowOrig = [0, 8000];   % window in which the IPSC response would lie (ms), original
ipscStartOrig = 1000;   % time of IPSC application (ms), original
ipscDur = 7000;         % duration of IPSC application (ms), for fitting

% Parameters for simulations
modelName = 'm3ha';   % 'm3ha' or 'soplata' or 'howard'
simMode = 'passive';    % 'passive' or 'active'
passiveOnly = true; %false;    % whether to include passive parameters only
closedLoop = false;     % whether to use the final state as the initial
                        %   condition for the next simulation
solverOrder = 0;        % uses the implicit Crank-Nicholson scheme
                        %   for multi-compartment models
timeToStabilize = 1000; % time for state variables to stabilize in ms

% Parameters for calculations and plotting
switch simMode
case 'passive'
    baseWindowOrig = [cprWindowOrig(1), cpDelayOrig];
    fitWindowOrig = [cpDelayOrig, cprWindowOrig(2)];
    xLimitsOrig = [cpDelayOrig - 10, cpDelayOrig + 40];
case 'active'
    baseWindowOrig = [ipscrWindowOrig(1), ipscStartOrig];
    fitWindowOrig = [ipscStartOrig, ipscrWindowOrig(2)];
    xLimitsOrig = [ipscStartOrig - 200, ipscStartOrig + 2000];
end

% Parameters for saving things
archiveScripts = false; % true;
archiveTopOnly = true;
saveModel = false; % true;
saveFigure = true;
createImportLog = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Generate full file paths
neuronParamsPath = fullfile(neuronParamsDir, neuronParamsFileName);

% Shift original times by timeToStabilize
cprWindow = cprWindowOrig + timeToStabilize;
ipscrWindow = ipscrWindowOrig + timeToStabilize;
cpDelay = cpDelayOrig + timeToStabilize;
fitWindow = fitWindowOrig + timeToStabilize;
baseWindow = baseWindowOrig + timeToStabilize;
xLimits = xLimitsOrig + timeToStabilize;

% Obtain the simulation end time in ms
switch simMode
case 'passive'
    timeEnd = cprWindow(2);
case 'active'
    timeEnd = ipscrWindow(2);
end

% Create file prefix
filePrefix = [create_time_stamp, '_', modelName, '_', simMode, '_', sweepName];

% Create output model file name
matFileName = [filePrefix, '_before_simulations.mat'];
figName = [filePrefix, '_comparison.png'];
if archiveTopOnly
    archiveName = [filePrefix, '_dependent_scripts_toponly.zip'];
else
    archiveName = [filePrefix, '_dependent_scripts.zip'];
end

% Create full paths
matFilePath = fullfile(outFolder, matFileName);
figPath = fullfile(outFolder, figName);

%% Import data to compare against
switch simMode
case 'passive'
    [sweepData, sweepInfo] = ...
            m3ha_import_raw_traces(sweepName, 'Directory', matFilesDir, ...
                    'Verbose', true, 'CreateLog', createImportLog, ...
                    'ToParsePulse', true, 'ToCorrectDcSteps', true, ...
                    'ToAverageByVhold', true, 'OutFolder', outFolder, ...
                    'TimeToPad', timeToStabilize, ...
                    'ResponseWindow', cprWindowOrig, ...
                    'StimStartWindow', cpStartWindowOrig);
case 'active'
    [sweepData, sweepInfo] = ...
            m3ha_import_raw_traces(sweepName, 'Directory', matFilesDir, ...
                    'Verbose', true, 'CreateLog', createImportLog, ...
                    'ToParsePulse', false, 'ToCorrectDcSteps', false, ...
                    'ToAverageByVhold', false, 'OutFolder', outFolder, ...
                    'TimeToPad', timeToStabilize, ...
                    'ResponseWindow', ipscrWindowOrig, ...
                    'StimStartWindow', ipscStartOrig);
end

% Extract from cell array
sweepData = sweepData{1};

% Extract needed data and information
[timeData, realVoltageData, realStimCopy] = extract_columns(sweepData, 1:3);
holdingPotential = sweepInfo.holdPotential;

% Update the sampling interval (esp. the number of rows)
siMs = compute_sampling_interval(timeData);

% Match the simulation time step with the sampling interval
timeStep = siMs;         % output time step in ms
simTimeStep = siMs;      % simulation time step in ms

%% Create the neuron
% Purge all pre-compiled binaries
xolotl.cleanup;

% Create a xolotl object based on a parameters file
%   Note: m3ha is a handle to the xolotl object
switch modelName
    case 'm3ha'
        m3ha = m3ha_xolotl_create_neuron(neuronParamsPath, ...
                                            'PassiveOnly', passiveOnly);
    case 'soplata'
        m3ha = xolotl_create_model_soplata;
    case 'howard'
        m3ha = xolotl_create_model_howard;
    otherwise
        error('modelName unrecognized!');
end

% Parse the xolotl object
parsedParams = parse_xolotl_object(m3ha);

% Get the number of compartments
nCompartments = parsedParams.nCompartments;

% Make sure holdingPotential has the correct dimensions
holdingPotential = match_dimensions(holdingPotential, [1, nCompartments]);

% Create NaN data for dendrites
nanDataDendrites = nan(size(realVoltageData, 1), nCompartments - 1);

% Put together as data to compare
dataToCompare = [nanDataDendrites, realVoltageData, realStimCopy];

% Set general simulation parameters
xolotl_set_simparams(m3ha, 'ClosedLoop', closedLoop, ...
                    'SolverOrder', solverOrder, 'Temperature', temperature, ...
                    'TimeStep', timeStep, 'SimTimeStep', simTimeStep);

% Set simulation parameters for the current pulse response
xolotl_set_simparams(m3ha, 'TimeEnd', timeEnd, ...
                        'InitialVoltage', holdingPotential);

% Add a current pulse protocol
%   Note: must do this after the correct t_end is set
switch simMode
case 'passive'
    xolotl_add_current_pulse(m3ha, 'Compartment', compToPatch, ...
                                'Delay', cpDelay, 'Duration', cpDuration, ...
                                'Amplitude', cpAmplitude);
    % xolotl_add_current_injection(m3ha, 'Compartment', compToPatch, ...
    %                             'CurrentVector', realStimCopy);
case 'active'
    xolotl_add_current_injection(m3ha, 'Compartment', compToPatch, ...
                                'CurrentVector', realStimCopy);
end

% Save the xolotl object before simulations
if saveModel
    save(matFilePath, 'm3ha');
end

%% Simulate and plot
% Simulate and use default plot
% m3ha.plot;

% Simulate and plot individual traces against data
m3haAfter = ...
    m3ha_xolotl_plot(m3ha, 'DataToCompare', dataToCompare, ...
                        'XLimits', xLimits, ...
                        'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
                        'FigTitle', figTitle, 'CompToPatch', compToPatch, ...
                        'TimeToStabilize', timeToStabilize, ...
                        'HoldingPotential', holdingPotential);

if saveFigure
    saveas(gcf, figPath);
end

%% Manipulate

% Set up the plots to be manipulated
m3ha.manipulate_plot_func = ...
    {@(x) m3ha_xolotl_plot(x, 'DataToCompare', dataToCompare, ...
                            'XLimits', xLimits, ...
                            'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
                            'FigTitle', figTitle, 'CompToPatch', compToPatch, ...
                            'TimeToStabilize', timeToStabilize)};

% Manipulate leak channel parameters
% m3ha.manipulate({'soma.len', 'soma.radius', 'soma.Leak.gbar', ...
%                   'dend1.len', 'dend1.radius', ...
%                   'dend2.len', 'dend2.radius'})
% m3ha.manipulate('*gbar')
% m3ha.manipulate('*Leak*')
m3ha.manipulate('*E')
% m3ha.manipulate({'*Leak*', '*length*'})
manip = gcf;
manip.Position(2) = manip.Position(2) - 200;


% Displays a list of properties
% properties(xolotl)

if archiveScripts
    archive_dependent_scripts(mfilename, 'OutFolder', outFolder, ...
                                'OutFileName', archiveName, ...
                                'TopOnly', archiveTopOnly);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

initialVoltageVc = -70; % initial voltage level for voltage clamp
timeEndVc = 1000;       % simulation end time in ms for voltage clamp

initialVoltageCpr = -60;% initial voltage level for current pulse response
m3ha = xolotl_set_simparams(m3ha, 'InitialVoltage', initialVoltageCpr, ...
                            'TimeEnd', timeEnd);

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
holdingPotential = compute_stats(vvecOrig, 'EndPoints', endPointsBase);

% Repeat the current pulse responses for each compartment for now
vvecsCpr = repmat(vvecCpr, [1, nCompartments]);

% Construct padding ones or zeros
vvecsToPadCpr = ones(nSamplesToPadForStabilization, 1) * holdingPotential;
ivecToPadCpr = zeros(nSamplesToPadForStabilization, 1);

realVoltageData = [vvecsToPadCpr; vvecsCpr];
realStimCopy = [ivecToPadCpr; ivecCpr];

dataToCompare = [realVoltageData, realVoltageData, ...
                realVoltageData, realStimCopy];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
