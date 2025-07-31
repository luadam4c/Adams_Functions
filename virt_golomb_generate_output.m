%% Runs Golomb 2022 Julia Code and Generates Output Figures
%
% File History:
% 2025-07-30 Created by Adam Lu
% 
% Requires:
%       cd/create_subplots.m
%       cd/plot_vertical_shade.m
%       cd/save_all_figtypes.m

% Adapted from notes from Matthew Bergosh: 
%   The basic structure of the model
% is shown in Fig. 2A of the paper. In the model, there are three
% populations of 100 neurons each:
%       vIRT retraction (I1) 
%       vIRT protraction (I2) 
%       facial motor neurons (FMN)
% The activity of the whiskers is then calculated based on activity in the
% FMN.
%   Constant current is applied to all three populations, and in the basic
% simulation I1 receives a randomly timed current input from the
% Pre-Botzinger complex (PB). The simulation calculates the voltages and
% spikes of all the neurons in response to their inputs, and iterates over
% time.
% File types:
%       irt_const.jl: stores parameters for each model
% Subfolders under dir_here:
%       pb_v_rast_intra: reflects the model shown in Fig. 2A
%       pb_v_rast_no_intra: I1-I1, I2-I2, I1-I2, I2-I1 all turned off
%
% irt.col Column Legend:
%       1 Time (ms)
%       2 Vav_pop[1]
%       3 sav_pop[1][1]
%       4 sav_pop[1][2]
%       5 Isyn_cont[1]
%       6 Varbar[1]
%       7 Varbar[8]
%       8 Vav_pop[2]
%       9 sav_pop[2][1]
%       10 sav_pop[2][2]
%       11 Isyn_cont[101]
%       12 Varbar[101]
%       13 Varbar[102]
%       14 Vav_pop[3]
%       15 sav_pop[3][1]
%       16 sav_pop[3][2]
%       17 Isyn_cont[201]
%       18 Varbar[201]
%   In brackets after the variable name is the population or neuron it's
% coming from. The first six columns are for the I1 pop, the second six for
% the I2, and the last five are for the FMN pop. 
%       Vpop is the average membrane potential across the population.
%       Sav_pop is unclear, some sort of synaptic transmission value. 
%           Based on synstr.synvar, not sure what that is.
%       Isyn seems like the total current coming into a population
%       Varb is some sort of voltage for a specific example neuron
%
% irt.sex Column Legend
%       1 Time (ms)
%       2 Active or not (binary)
%
% NEW SIMULATIONS
%   When you want to create a new simulation, create a folder with a
% modified copy of irt_const.jl. Then, in the irt.jl file under "prog", 
% change line269 to reflect the title of the new folder, and run the code 
% using a Julia REPL. This will execute the simulation and save the output 
% files in the new folder.

%% Hard-coded parameters
packagesNeeded = ["Formatting", "SmoothingSplines", ...
                    "Polynomials", "DSP", ...
                    "NLsolve", "StatsBase"];
dirData = 'datfig';
dirSubData = 'dir_here';
dirFig = 'figMCode';
dirIntra = 'pb_v_rast_intra';
dirNoIntra = 'pb_v_rast_no_intra';
scriptInput = fullfile('genfig', 'gen_input_files.jl');
scriptIrt = fullfile('prog', 'irt.jl');
fileToCheck = 'run_here.com';
fileVolt = 'irt.col';   % voltages of populations and example neurons 
filePBotC = 'irt.sex';  % pre-Botzinger Complex conductance trace
fileFig6A = 'fig6A';

% Model parameters
gInput = 0.5;           % Input conductance from pre-BotC to I1 (mS/cm2)

% Plotting parameters
timeLim = [3.5, 6];    % Time (s)
gLim = [-0.05, 0.55];   % Conductance (mS/cm2)
vLim = [-75, 60];       % Membrane potential (mV)

% Save parent directory
pathParent = pwd;

%% PREPARATION
% Check if julia is installed
% Install from https://julialang.org/downloads/
juliaInstalled = system(['where ', 'julia']) == 0;
if ~juliaInstalled
    disp('Install julia from https://julialang.org/downloads/');
    return
end

% Check if julia packages are installed, and provide instructions if not:
% TODO
% Packages needed: 
%   Formatting SmoothingSplines Polynomials DSP NLsolve StatsBase
for iPackage = 1:length(packagesNeeded)
    msg = sprintf('Within julia please enter ''import Pkg; Pkg.add("%s")''', ...
                    packagesNeeded(iPackage));
    disp(msg)
end

% Check if directory contains data directory
if ~isfolder(dirData)
    fprintf("Please run this code with the '%s' folder in current directory\n", dirData);
    return
end

%% Run Simulations
% Change directory to data directory
cd(dirData);

% Generate input files
pathToCheck = fullfile(pathParent, dirData, fileToCheck);
if exist(pathToCheck, 'file') ~= 2
    command0 = sprintf('julia %s', fullfile(pathParent, scriptInput));
    system(command0);
else
    fprintf('Input files already generated. %s already exists!\n', pathToCheck);    
end

% Change directory to data subdirectory
cd(dirSubData);
pathSubData = fullfile(pathParent, dirData, dirSubData);

% Generate data for figures 6A
pathVoltIntra = fullfile(pathSubData, dirIntra, fileVolt);
if exist(pathVoltIntra, 'file') ~= 2
    command1 = sprintf('julia %s %s', fullfile(pathParent, scriptIrt), dirIntra);
    system(command1);
else
    fprintf('%s already exists!\n', pathVoltIntra);
end

% Generate data for figures ???
pathVoltNoIntra = fullfile(pathSubData, dirNoIntra, fileVolt);
if exist(pathVoltNoIntra, 'file') ~= 2
    command1 = sprintf('julia %s %s', fullfile(pathParent, scriptIrt), dirNoIntra);
    system(command1);
else
    fprintf('%s already exists!\n', pathVoltNoIntra);
end

%% Create Figure 6A
% Create paths
pathPBIntra = fullfile(pathSubData, dirIntra, filePBotC);

% Read Data for figures with intrapopulation connections
if ~exist('voltIntra', 'var')
    voltIntra = load(pathVoltIntra);
end
if ~exist('pbIntra', 'var')
    pbIntra = load(pathPBIntra);
end

% Read time vectors
tIntra = voltIntra(:, 1)/1000;
tPbIntra = pbIntra(:, 1)/1000;

% Read input conductance data;
gPbIntra = gInput * pbIntra(:, 2);

% Extract PreBotz On and Off times
isOn = gPbIntra ~= 0;
isOff = gPbIntra == 0;
isPrevOn = [0; isOn(1:end-1)];
isNextOn = [isOn(2:end); 0];
isStart = isOn & isNextOn;
isEnd = isOn & isPrevOn;
timesStart = tPbIntra(isStart);
timesEnd = tPbIntra(isEnd);
nStims = sum(isStart);
timesStim = cell(nStims, 1);
for iStim = 1:nStims
    timesStim{iStim} = [timesStart(iStim), timesEnd(iStim)];
end

% Read individual neuron data
vIntraI1Cell1 = voltIntra(:, 6);
vIntraI1Cell2 = voltIntra(:, 7);
vIntraI2Cell1 = voltIntra(:, 12);
vIntraI2Cell2 = voltIntra(:, 13);
vIntraFmnCell1 = voltIntra(:, 18);

% Read population average data
vVecI1 = voltIntra(:, 2);
vVecI2 = voltIntra(:, 8);
vVecFMN = voltIntra(:, 14);

% Create subplots
[fig6A, ax6A] = create_subplots(4, 1, 'FigExpansion', [0.3, 1]);

% Plot input conductance data
axes(ax6A(1));
plot(tPbIntra, gPbIntra);
plot_vertical_shade(timesStim);
xlim(timeLim);
ylim(gLim);
ylabel('g_{rB} (mS/cm^2)');

% Plot voltage of I1Cell1
axes(ax6A(2));
plot(tIntra, vIntraI1Cell1);
plot_vertical_shade(timesStim);
xlim(timeLim);
ylim(vLim);
ylabel('vIRt_{ret} (mV)');

% Plot voltage of I2Cell1
axes(ax6A(3));
plot(tIntra, vIntraI2Cell1);
plot_vertical_shade(timesStim);
xlim(timeLim);
ylim(vLim);
ylabel('vIRt_{pro} (mV)');

% Plot voltage of FmnCell1
axes(ax6A(4));
plot(tIntra, vIntraFmnCell1);
plot_vertical_shade(timesStim);
xlim(timeLim);
ylim(vLim);
ylabel('vFMN (mV)');

figure(fig6A);
xlabel('Time (s)');

% Save to file
pathFig6A = fullfile(pathParent, dirFig, fileFig6A);
save_all_figtypes(fig6A, pathFig6A);

%% Create Figure 6B
% TODO

% Read Data for figures without intrapopulation connections
% voltNoIntra = load(pathVoltNoIntra);  % Load as an array