function m3ha_network_launch(nCells, useHH, candidateIDs)
%% Launches NEURON with simulation commands and plot output figures
%
% Requires:
%       cd/argfun.m
%       cd/check_dir.m
%       cd/compile_mod_files.m
%       cd/construct_fullpath.m
%       cd/create_labels_from_numbers.m
%       cd/create_looped_params.m
%       cd/create_time_stamp.m
%       cd/find_in_strings.m
%       cd/load_params.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_network_change_params.m
%       cd/m3ha_network_show_net.m
%       cd/m3ha_network_raster_plot.m
%       cd/m3ha_network_single_neuron.m
%       cd/m3ha_network_define_actmode.m
%       cd/run_neuron.m
%       cd/save_params.m
%       /media/adamX/m3ha/network_model/m3ha_run_1cell.hoc
%       /media/adamX/m3ha/network_model/m3ha_net.hoc
%       /media/adamX/m3ha/network_model/RE.tem
%       /media/adamX/m3ha/network_model/TC3.tem
%       /media/adamX/m3ha/network_model/ipulse2.mod
%       /media/adamX/m3ha/network_model/gabaA_Cl.mod
%       /media/adamX/m3ha/network_model/gabaa.mod (potentially)

% File History:
% 2017-10-23 Modified from /RTCl/neuronlaunch112.m
% 2017-10-30 Added pCond, gIncr, etc.
% 2017-10-31 Removed REuseca and added useHH
% 2017-11-03 Added bicucullineFlag & heteroTCFlag
% 2017-11-04 gIncr now scales TCgabaa as well
% 2017-11-06 Moved code to define_actmode.m
% 2017-11-07 Added candidateIDs
% 2017-11-07 Added actMode = 8~10
% 2017-11-08 Added templateLabel to fileLabel
% 2017-11-08 Added seedNumber
% 2018-02-28 Don't use HH
% 2018-03-29 Fixed ordering of useHH and REnsegs
% 2018-04-17 Plot inSlopeWatching if useHH is 0
% 2018-04-24 Moved pfiles sources from ../optimizer4gabab/pfiles 
%               to just pfiles
% 2018-04-26 Added the case for nCells == 2
% 2018-04-26 Made candidateIDs an argument
% 2018-04-26 Made nCells and useHH arguments
% 2019-10-31 Now uses m3ha_net.hoc
% 2019-10-31 Now uses run_neuron.m
% 2019-11-04 Updated for m3ha_network_loop_launch_20191104.m
% 2019-11-08 Fixed RERErad, TCRErad, RETCrad
% 2019-11-18 Added modifications for running on Windows
% TODO: Plot gAMPA and gGABA instead of the i's for synaptic event monitoring
% TODO: Make the network circular to lose edge effects
% TODO: Perform simulations to generate a linear model
% TODO: Update specs for m3ha_network_raster_plot.m
%

%% Arguments defined here temporarily
% nCells = 1;
% nCells = 2;
% nCells = 20;
% nCells = 100;
% useHH = true;           % whether to use HH channels
% candidateIDs = [25, 36, 27];
% candidateIDs = 25;
% candidateIDs = 22;

%% Experiment Name
experimentName = 'm3ha';
% experimentSuffix = 'stimstart_3000';
experimentSuffix = ['ncells_', num2str(nCells), '_useHH_', num2str(useHH), ...
            '_templateIDs_', strjoin(strsplit(num2str(candidateIDs)), ',')];

%% Hard-coded parameters
homeDirName = 'network_model';
paramsDirName = fullfile('best_params', ...
                'bestparams_20171213_singleneuronfitting16_Rivanna');
% paramsDirName = fullfile('best_params', ...
%                 'bestparams_20180424_singleneuronfitting21_Rivanna');

%% Flags
debugFlag = false;              % whether to do a very short simulation
seedNumber = 1;                 % number to seed random number generator
simNumbers = []; %5;            % run only these simulation numbers
onLargeMemFlag = false;         % whether to run on large memory nodes
onHpcFlag = false;              % whether on high-performance computing server
saveAllVariablesFlag = false;   % whether to save variables as a .mat file
saveStdOutFlag = false;         % whether to always save standard outputs
bicucullineFlag = true;         % whether GABA-A conductances are removed
loopMode = 'grid'; %cross;      % how to loop through parameters: 
                                %   'cross' - Loop through each parameter 
                                %               while fixing others
                                %   'grid'  - Loop through all possible 
                                %               combinations of parameters

% Decide on what to save and plot
if nCells == 1 || nCells == 2
    savePlotMode = 'spikes&special';
elseif nCells == 20 || nCells == 100
    savePlotMode = 'spikes';    
else
    error('nCells = %d is not implemented yet!', nCells);
end
savePlotMode = 'spikes';

%% Simulation modes
simMode = 1;    % 1 - full simulation
                % 2 - short simulation
                % 3 - medium simulation

%% Activation modes
% 1 - Activate a single RE cell by injecting a train of current pulses
% 2 - Activate every (RERErad + 1)th RE cell by injecting trains of 
%       current pulses
% 3 - Activate 3 RE cells (RERErad + 1) cells apart by injecting 
%       trains of current pulses
% 4 - Activate a single RE cell by changing the membrane potential 
%       instantaneously
% 5 - Activate RE cells with a Gaussian likelihood by changing 
%       the membrane potential instantaneously
% 6 - Activate every 3rd RE cell by injecting trains of current pulses
% 7 - Activate all RE cells by injecting trains of current pulses
% 8 - Activate 3 RE cells RETCrad cells apart by injecting trains of 
%       current pulses
% 9 - Activate 10 center RE cells by injecting current pulses
% 10 - Activate 20 center RE cells by injecting current pulses

% actMode = 10;   
actMode = 1;   

% Decide on template TC neurons to use
%% Template TC neurons;
templateNames = {'D091710'; 'E091710'; 'B091810'; 'D091810'; ...
                'E091810'; 'F091810'; 'A092110'; 'C092110'; ...
                'B092710'; 'C092710'; 'E092710'; 'A092810'; ...
                'C092810'; 'K092810'; 'A092910'; 'C092910'; ...
                'D092910'; 'E092910'; 'B100110'; 'E100110'; ...
                'A100810'; 'B100810'; 'D100810'; 'A101210'; ...
                'C101210'; 'D101210'; 'E101210'; 'F101210'; ...
                'I101210'; 'M101210'; 'B101310'; 'D101310'; ...
                'E101310'; 'F101310'; 'G101310'; 'H101310'};
%candidateIDs = 21;
%candidateIDs = 3;
%candidateIDs = 20;
%candidateIDs = 21;
%candidateIDs = 22;
%candidateIDs = 26;
%candidateIDs = 27;
%candidateIDs = 33;
%candidateIDs = 25;
%candidateIDs = [22; 27];
%candidateIDs = [22; 33];
%candidateIDs = [20; 22; 33];
%candidateIDs = [20; 22; 27; 33];
%candidateIDs = [3; 20; 22; 26; 27; 33];
%candidateIDs = [3; 20; 22; 25; 26; 27; 33];
%candidateIDs = [3; 20; 21; 22; 25; 26; 27; 33];
%candidateIDs = 1:36;
%candidateIDs = [22; 27; 33];      % 'B100810', 'E101210', 'E101310'

% In ascending order of total error in singleneuronfitting16_Rivanna:
%candidateIDs = 25;                 % 'C101210'
%candidateIDs = 36;                 % 'H101310'
%candidateIDs = 27;                 % 'E101210'
%candidateIDs = [25; 36; 27];       % 'C101210', 'H101310', 'E101210'

%candidateIDs = 20;                 % 'E100110'
%candidateIDs = 22;                 % 'B100810'
%candidateIDs = 30;                 % 'M101210'
%candidateIDs = 34;                 % 'F101310'
%candidateIDs = 10;                 % 'C092710'
%candidateIDs = 21;                 % 'A100810'
%candidateIDs = 3;                  % 'B091810'

%% hoc file names
switch nCells
    case 100
        hocFileName = 'm3ha_run_100cells.hoc';
    case 20
        hocFileName = 'm3ha_run_20cells.hoc';
    case 2
        hocFileName = 'm3ha_run_2cells.hoc';
    case 1
        hocFileName = 'm3ha_run_1cell.hoc';
    otherwise
        error('nCells == %d is not implemented yet!', nCells);
end

%% For parpool
if onLargeMemFlag || debugFlag || numel(simNumbers) == 1 || ...
        simMode == 2 || nCells <= 2
    % No need renew parpool each batch if memory is not an issue
    renewParpoolFlagNeuron = 0;    % whether to renew parpool every batch to release memory
    maxNumWorkersNeuron = 20;      % maximum number of workers for running NEURON 
    renewParpoolFlagPlots = 0;     % whether to renew parpool every batch to release memory
    maxNumWorkersPlots = 20;       % maximum number of workers for plotting things
else
    switch savePlotMode
    case 'all'          % saving and plotting everything
        %% For parpool
        renewParpoolFlagNeuron = 1;% whether to renew parpool every batch to release memory
        maxNumWorkersNeuron = 12;  % maximum number of workers for running NEURON 
        renewParpoolFlagPlots = 1; % whether to renew parpool every batch to release memory
        maxNumWorkersPlots = 12;   % maximum number of workers for plotting things
    case 'curves'           % saving spikes and plotting curves/maps only
        renewParpoolFlagNeuron = 0;% whether to renew parpool every batch to release memory
        maxNumWorkersNeuron = 20;  % maximum number of workers for running NEURON 
        renewParpoolFlagPlots = 0; % whether to renew parpool every batch to release memory
        maxNumWorkersPlots = 20;   % maximum number of workers for plotting things
    case 'spikes'           % saving spikes and plotting raster plots and curves/maps only
        renewParpoolFlagNeuron = 0;% whether to renew parpool every batch to release memory
        maxNumWorkersNeuron = 20;  % maximum number of workers for running NEURON 
        renewParpoolFlagPlots = 0; % whether to renew parpool every batch to release memory
        maxNumWorkersPlots = 20;   % maximum number of workers for plotting things
    case 'spikes&special'   % saving spikes and special neuron traces only
        renewParpoolFlagNeuron = 0;% whether to renew parpool every batch to release memory
        maxNumWorkersNeuron = 12;  % maximum number of workers for running NEURON 
        renewParpoolFlagPlots = 1; % whether to renew parpool every batch to release memory
        maxNumWorkersPlots = 12;   % maximum number of workers for plotting things
    end
end

switch savePlotMode
case 'all'
    %% Save flags
    saveNetwork = 1;        % whether to save network topology
    saveSpikes = 1;         % whether to save spike data
    saveSomaVoltage = 1;    % whether to save all voltage data
    saveSomaCli = 1;        % whether to save all chloride concentration data
    saveSpecial = 1;        % whether to save special neuron data

    %% Plot flags
    plotNetwork = 1;        % whether to plot network topology
    plotSpikes = 1;         % whether to plot spike data
    plotTuning = 1;         % whether to plot tuning curves
    plotSingleNeuronData = 1;% whether to plot single neuron data
case 'curves'
    %% Save flags
    saveNetwork = 0;        % whether to save network topology
    saveSpikes = 1;         % whether to save spike data
    saveSomaVoltage = 0;    % whether to save all voltage data
    saveSomaCli = 0;        % whether to save all chloride concentration data
    saveSpecial = 0;        % whether to save special neuron data

    %% Plot flags
    plotNetwork = 0;        % whether to plot network topology
    plotSpikes = 0;         % whether to plot spike data
    plotTuning = 1;         % whether to plot tuning curves
    plotSingleNeuronData = 0;% whether to plot single neuron data
case 'spikes'
    %% Save flags
    saveNetwork = 0;        % whether to save network topology
    saveSpikes = 1;         % whether to save spike data
    saveSomaVoltage = 0;    % whether to save all voltage data
    saveSomaCli = 0;        % whether to save all chloride concentration data
    saveSpecial = 0;        % whether to save special neuron data

    %% Plot flags
    plotNetwork = 0;        % whether to plot network topology
    plotSpikes = 1;         % whether to plot spike data
    plotTuning = 1;         % whether to plot tuning curves
    plotSingleNeuronData = 0;% whether to plot single neuron data
case 'spikes&special'
    %% Save flags
    saveNetwork = 0;        % whether to save network topology
    saveSpikes = 1;         % whether to save spike data
    saveSomaVoltage = 0;    % whether to save all voltage data
    saveSomaCli = 0;        % whether to save all chloride concentration data
    saveSpecial = 1;        % whether to save special neuron data

    %% Plot flags
    plotNetwork = 0;        % whether to plot network topology
    plotSpikes = 1;         % whether to plot spike data
    plotTuning = 1;         % whether to plot tuning curves
    plotSingleNeuronData = 1;% whether to plot single neuron data
end

% Code not fixed yet
plotTuning = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters to loop through
%{
pNames  = {'REdiam'};       % names of parameters to loop through
pLabels = {'REdiam (um)'};  % labels of parameters to loop through
pMin    = [2];              % minimum values of parameters to loop through
pMax    = [20];             % maximum values of parameters to loop through
pInc    = [2];              % increments of parameters to loop through
pIsLog  = [0];              % whether increments of parameters is in log
%}
%{
pNames  = {'REgabaGmax'};      % names of parameters to loop through
pLabels = {'REgabaGmax (uS)'}; % labels of parameters to loop through
pMin    = [0.0025];         % minimum values of parameters to loop through
pMax    = [0.045];          % maximum values of parameters to loop through
pInc    = [0.0025];         % increments of parameters to loop through
pIsLog  = [0];              % whether increments of parameters is in log
%}
%{
pNames  = {'stimFreq'};    % names of parameters to loop through
pLabels = {'Stimulation Frequency (Hz)'};% labels of parameters to loop through
pMin    = [1];              % minimum values of parameters to loop through
pMax    = [128];            % maximum values of parameters to loop through
pInc    = [2^(1/2)];        % increments of parameters to loop through
pIsLog  = [1];              % whether increments of parameters is in log
%}
%{
pNames  = {'REtauKCC2'};    % names of parameters to loop through
pLabels = {'Time constant of KCC2 (s)'};    % labels of parameters to loop through
pMin    = [4];              % minimum values of parameters to loop through
pMax    = [64];             % maximum values of parameters to loop through
pInc    = [2^(1/4)];        % increments of parameters to loop through
pIsLog  = [1];              % whether increments of parameters is in log
%}
%{
pNames  = {'REdiam', 'REgabaGmax', 'stimFreq', 'REtauKCC2'};    % names of parameters to loop through
pLabels = {'REdiam (um)', 'REgabaGmax (uS)', 'Stimulation Frequency (Hz)', 'Time constant of KCC2 (s)'};    % labels of parameters to loop through
pMin    = [4, 0.0025, 0.125, 0.25]; % minimum values of parameters to loop through
pMax    = [15, 0.045, 256, 64];     % maximum values of parameters to loop through
pInc    = [1, 0.0025, 2, sqrt(2)];  % increments of parameters to loop through
pIsLog  = [0, 0, 1, 1];             % whether increments of parameters is in log
%}
%{
pNames  = {'stimFreq', 'REtauKCC2'};   % names of parameters to loop through
pLabels = {'Stimulation Frequency (Hz)', 'Time constant of KCC2 (s)'};    % labels of parameters to loop through
pMin    = [16, 32*2^(1/4)];         % minimum values of parameters to loop through
pMax    = [19, 64];                 % maximum values of parameters to loop through
pInc    = [0.1, 2^(1/64)];          % increments of parameters to loop through
pIsLog  = [0, 1];                   % whether increments of parameters is in log
%}
%{
pNames  = {'REdiam', 'REgabaGmax'};    % names of parameters to loop through
pLabels = {'REdiam (um)', 'REgabaGmax (uS)'};    % labels of parameters to loop through
pMin    = [8, 0.1];                 % minimum values of parameters to loop through
pMax    = [12, 0.5];                % maximum values of parameters to loop through
pInc    = [0.5, 0.05];              % increments of parameters to loop through
pIsLog  = [0, 0];                   % whether increments of parameters is in log
%}

%% Parameters for various GABA-B conductance profiles
% pCond = 1;      % Pharmacological condition
                %   1 - Control; 2 - GAT 1 Block; 3 - GAT 3 Block; 4 - Dual Block
% gIncr = 100;    % GABA-B conductance amplitude scaling (%)

%{
pCond = 1;
gIncr = 60;
pNames  = {'pCond', 'gIncr'};       % names of parameters to loop through
pLabels = {'Pharm Condition', 'gGABAB amp scaling (%)'};  % labels of parameters to loop through
pMin    = [1, 30];                  % minimum values of parameters to loop through
pMax    = [4, 90];                  % maximum values of parameters to loop through
pInc    = [1, 30];                  % increments of parameters to loop through
pIsLog  = [0, 0];                   % whether increments of parameters is in log
%}

%{
pCond = 1;
gIncr = 15;
pNames  = {'pCond', 'gIncr'};       % names of parameters to loop through
pLabels = {'Pharm Condition', 'gGABAB amp scaling (%)'};  % labels of parameters to loop through
pMin    = [1, 7.5];                 % minimum values of parameters to loop through
pMax    = [4, 22.5];                % maximum values of parameters to loop through
pInc    = [1, 7.5];                 % increments of parameters to loop through
pIsLog  = [0, 0];                   % whether increments of parameters is in log
%}

pCond = 1;
gIncr = 30;
pNames  = {'pCond', 'gIncr'};       % names of parameters to loop through
pLabels = {'Pharm Condition', 'gGABAB amp scaling (%)'};  % labels of parameters to loop through
pMin    = [1, 5];                   % minimum values of parameters to loop through
pMax    = [4, 60];                  % maximum values of parameters to loop through
pInc    = [1, 5];                   % increments of parameters to loop through
pIsLog  = [0, 0];                   % whether increments of parameters is in log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global parameters to be defined at the start of NEURON
%   Will call either m3ha_run5.hoc or m3ha_run5_small.hoc or m3ha_run12_small.hoc
% nCells = 2; %100; %20;

%% Global parameters to be defined at the start of NEURON, to be consistent with m3ha_run5.hoc
celsius = 33;   % SET IN m3ha_run5.hoc: temperature of experiment (celsius) % TODO: What did Mark use?
nSpecial = 2;        % SET IN m3ha_run5.hoc: number of special RE or TC neurons
if nCells == 100 || nCells == 20
    RERErad = 8;    % SET IN m3ha_run5.hoc: radius of intra-RE connections
                    %   Sohal & Huguenard 2004 used 8
                    %   RTCl & Sohal & Huguenard 2003 used 4
elseif nCells == 2
    RERErad = 1;
elseif nCells == 1
    RERErad = 0;
else
    error('nCells unrecognized!');
end

%% Network parameters
if nCells == 100 || nCells == 20
    TCRErad = 4; %1; %2;    % radius of TC-RE connections
                    %   Sohal & Huguenard 2004 used 2
    RETCrad = 8; %4; % radius of RE-TC connections
                    %   Sohal & Huguenard 2004 used 4
elseif nCells == 2
    TCRErad = 1;
    RETCrad = 1;
elseif nCells == 1
    TCRErad = 0;
    RETCrad = 0;
else
    error('nCells unrecognized!');
end
spThr = 0;     % action potential threshold (mV)
synDel = 1;    % synaptic delay (ms)
synWeight = 1;      % synaptic weight (fraction of channels activated)
                %     for simplicity assume channels are always activated and that 
                %     channels have linearly additive effects


%% RE neuron parameters
REnsegs = 1;    % number of segments in an RE cell (1, 3 or 9)
                %     if REnsegs >= 3, the GABAA synapses will be distributed on each side
REcldnum = 2;   % which cld mechanism to use (0 or 1 or 2) in RE cells
REconsyn = 0;   % whether to concentrate synapses (0 or 1) in RE cells
                %     for REnsegs = 1 or 3, we have soma, soma_flank[0], soma_flank[1]
REtauKCC2 = 32; % Cl- removal time constant (s) in RE cells
REepas = -70;   % leak reversal potential (mV) of RE cells, Peter's value
                %    Sohal & Huguenard 2003 & 2004 used -77 mV
                %    Jedlicka et al 2011 used -60 mV
REdiam = 10;    % diameter (um) of an RE cell, Peter's value
REgpasLB = 4.5e-5;  % lower bound for passive leak conductance (S/cm^2) in RE cells, Sohal & Huguenard 2003
REgpasUB = 5.5e-5;  % upper bound for passive leak conductance (S/cm^2) in RE cells, Sohal & Huguenard 2003
                %     Jedlicka et al 2011 used 2e-4 S/cm^2 %%% What should we use?

%% Synapse parameters
if bicucullineFlag
    REgabaGmax = 0;
else
    REgabaGmax = 0.005; %0;     % maximal conductance (uS) of the GABA-A receptor on RE cells
                            %   Sohal & Huguenard 2003 varied between 5~12.5 nS
                            %   RTCl used 0.02 nS
end
%REampaGmax = 0.007 * (2*TCRErad + 1);
REampaGmax = 0.007;
                            % maximal conductance (uS) of the AMPA receptor on RE cells
                            %   Deleuze & Huguenard 2016 has about 7 nS per synapse 
                            %       (a minimal stimulation protocol was used)
                            %   Sohal & Huguenard 2004 used 0.05 uS
TCgabaaErev = -80;          % reversal potential (mV) of the GABA-A receptor on TC cells
                            %   Sohal & Huguenard 2004 used -85 mV; Traub 2005 used -81 mV
                            %   Peter's measurement gave -80 mV                           
if bicucullineFlag
    TCgabaaGmax = 0;
else
    TCgabaaGmax = 0.00896;  % maximal conductance (uS) of the GABA-A receptor on TC cells 
                            %   Sohal & Huguenard 2004 used 0.1 uS
                            %   The maximal GABA-B conductance is 0.00448 uS
                            %   Based on Huguenard & Prince, 1994, 
                            %       the maximal GABA-A conductance is about 5 times that of GABA-B away from Erev
                            %       and about 2 time that of GABA-B close to Erev
                            %   This must be consistent with /home/Matlab/Adams_Functions/m3ha_network_update_dependent_params.m
end
TCgababErev = -100; %-105; %-115        % reversal potential (mV) of the GABA-B receptor on TC cells
                            %   Christine used -115 mV in dynamic clamp experiments
                            %   Huguenard & Prince 1994 has -105 mV
                            %   ek is -100 mV
TCgababAmp = 0.016;         % conductance amplitude (uS) of the GABA-B receptor on TC cells
TCgababTrise = 52;          % rising phase time constant (ms) of the GABA-B receptor on TC cells
TCgababTfallFast = 90.1;    % fast decay time constant (ms) of the GABA-B receptor on TC cells
TCgababTfallSlow = 1073.2;  % slow decay time constant (ms) of the GABA-B receptor on TC cells
TCgababW = 0.952;           % weight (1) of the fast decay of the GABA-B receptor on TC cells

%% Initial ion concentrations
cai0 = 2.4e-4;  % initial intracellular [Ca++] (mM), Destexhe et al
cao0 = 2;       % initial extracellular [Ca++] (mM), Peter's value
% cli0 = 5;     % initial intracellular [Cl-] (mM), 
                %   corresponding to eGABA = -71 mV
                %            TODO: Ulrich
cli0 = 8;       % initial intracellular [Cl-] (mM),
                %   corresponding to eGABA = -61 mV
                %            Jedlicka et al 2011 (agrees with Peter's data)
% cli0 = 11;    % initial intracellular [Cl-] (mM),
                %   corresponding to eGABA = -54 mV
                %            TODO: Peter's data with ChABC added
% cli0 = 17;    % initial intracellular [Cl-] (mM),
                %   corresponding to eGABA = -45 mV
                %            TODO: Sun
clo0 = 130.5;   % initial extracellular [Cl-] (mM), Peter's value (Jedlicka et al 2011 used 133.5 mM)

%% Activation parameters for actMode == 1~3, 6~9
actCellID = floor(nCells/2);    % ID # of central neuron to activate
stimStart = 3000; %300;        % stimulation delay (ms)
stimDur = 40;                  % stimulation duration (ms)
stimFreq = 0.1;                % stimulation frequency (Hz),
                                %   must be less than 1000/cpDur
cpDur = 40;                    % current pulse duration (ms)

% The following must be consistent with m3ha_network_update_dependent_params.m
cpAmp = 0.2*(REdiam/10)^2;     % current pulse amplitude (nA),
                                %   must be proportional to square of diameter 
cpPer = floor(1000/stimFreq); % current pulse period (ms),
                                %   i.e. interval between pulse onsets
cpNum = ceil(stimDur/cpPer); % number of current pulses

%% Activation parameters for actMode == 4
actCellV = 0;                   % voltage (mV) to set activated neuron to

%% Activation parameters for actMode == 5
actWidth = 50;                  % width of Gaussian distribution for
                                %   randomly activating cells
actMaxP = 0.5;                  % maximum likelihood of activation at center

%% Simulation parameters
if simMode == 1
    tStop = 30000;              % total time of simulation (ms)
elseif simMode == 2
    tStop = 4000; %1000;        % total time of simulation (ms)
elseif simMode == 3
    tStop = 7000; %4000;        % total time of simulation (ms)
end
dt = 0.1;                       % time step of integration (ms)

%% Recording parameters
if simMode == 1
    tStart = 0;                 % time to start plotting (ms)
elseif simMode == 2
    tStart = 0;                 % time to start plotting (ms)
elseif simMode == 3
    tStart = 0;                 % time to start plotting (ms)
end
REsp1cellID = actCellID;        % ID # of 1st special RE neuron to record
if nCells > 1
    REsp2cellID = actCellID - 1;    % ID # of 2nd special RE neuron to record
else
    REsp2cellID = actCellID;    % ID # of 2nd special RE neuron to record
end
TCsp1cellID = actCellID;        % ID # of 1st special TC neuron to record
if nCells > 1
    TCsp2cellID = actCellID - 1;    % ID # of 2nd special TC neuron to record
else
    TCsp2cellID = actCellID;    % ID # of 2nd special TC neuron to record
end

%% Set ID #s of neurons to plot
act = actCellID;            % ID # of the activated neuron
if nCells == 100
    actLeft1 = actCellID - 1;  % ID # of the neuron one below the activated neuron
    actLeft2 = actCellID - 10; % ID # of the neuron 10 below the activated neuron
    far = actCellID - 20;       % ID # of a far away neuron
elseif nCells == 20
    actLeft1 = actCellID - 1;  % ID # of the neuron one below the activated neuron
    actLeft2 = actCellID - 2;  % ID # of the neuron 2 below the activated neuron
    far = actCellID - 10;       % ID # of a far away neuron
elseif nCells == 2
    actLeft1 = actCellID - 1;  % ID # of the neuron one below the activated neuron
    actLeft2 = actCellID;      % Repeat for compatibility
    far = actCellID - 1;        % Repeat for compatibility
elseif nCells == 1
    actLeft1 = actCellID;       % Repeat for compatibility
    actLeft2 = actCellID;       % Repeat for compatibility
    far = actCellID;            % Repeat for compatibility
else
    error('nCells unrecognized!');
end

%% Arguments for plotting (not logged in sim_params)
propertiesToPlot = 1:8;         % property #s of special neuron to record to be plotted (maximum range: 1~8, must be consistent with net.hoc)
%[1, 5, 6, 8]; 
cellsToPlot = [act, actLeft1, actLeft2, far]; % ID #s for neurons whose voltage is to be plotted

%% Set output file names; must have only one '.' (not logged in sim_params)
simparamsF = 'sim_params.csv';  % file with simulation parameters
scmdsF = 'sim_commands.txt';    % file with simulation commands
soutF = 'sim_output.txt';       % file with simulation standard outputs
sREREsynF = 'RERE.syn';         % file with RE-RE synaptic connections
sTCREsynF = 'TCRE.syn';         % file with TC-RE synaptic connections
sRETCsynF = 'RETC.syn';         % file with RE-TC synaptic connections
sREspikeF = 'RE.spi';           % file with RE spike train output
sTCspikeF = 'TC.spi';           % file with TC spike train output
sREvF = 'RE.singv';             % file with RE single neuron voltage traces
sTCvF = 'TC.singv';             % file with TC single neuron voltage traces
sREcliF = 'RE.singcli';         % file with RE single neuron chloride concentration traces
sREsp1F = ['RE[', num2str(REsp1cellID), '].singsp'];    % file with RE special neuron #1 other traces
sREsp2F = ['RE[', num2str(REsp2cellID), '].singsp'];    % file with RE special neuron #2 other traces
sTCsp1F = ['TC[', num2str(TCsp1cellID), '].singsp'];    % file with TC special neuron #1 other traces
sTCsp2F = ['TC[', num2str(TCsp2cellID), '].singsp'];    % file with TC special neuron #2 other traces
sREleakF = 'REleak.csv';        % file with RE neuron leak properties
sREparamsF = 'REparams.csv';    % file with RE neuron parameters
sTCparamsF = 'TCparams.csv';    % file with TC neuron parameters

%% For debug mode
if debugFlag
    tStart = 0;
    tStop = 2000;

    % Minimize number of points
    for p = 1:length(pNames)
        if pIsLog(p)
            pInc(p) = pMax(p)/pMin(p);
        else
            pInc(p) = pMax(p) - pMin(p);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of template candidates
nCandidates = length(candidateIDs);

% Set a flag for whether heterogeneity is introduced
if nCandidates > 1
    heteroTCFlag = 1;
else
    heteroTCFlag = 0;
end

%% Set folders for reading and saving files
% Find parent and home directory
parentDirectory = m3ha_locate_homedir;
homeDirectory = fullfile(parentDirectory, homeDirName);

% Compile or re-compile .mod files in the home directory
compile_mod_files(homeDirectory);

% Make directory to save all data
%   Note: Use current date & time in the format: YYYYMMDDThhmm
timeStamp = create_time_stamp('FormatOut', 'yyyymmddTHHMM');
if nCandidates > 1
    templateLabel = ['hetero', num2str(nCandidates)];
else
    templateLabel = templateNames{candidateIDs(1)};
end
fileLabel = [timeStamp, '_', templateLabel, '_', experimentSuffix];
outFolder = fullfile(homeDirectory, fileLabel);     
check_dir(outFolder);

%% Construct looped parameters
[pchnames, pchvalues, nSims, nump, pvalues, nperp] = ...
    create_looped_params (loopMode, pNames, pLabels, pIsLog, pMin, pMax, pInc, ...
            'OutFolder', outFolder, 'FileLabel', fileLabel, ...
            'NCells', nCells, 'ActMode', actMode);

% Force as a cell array
if isnumeric(pchvalues)
    pchvaluesCell = num2cell(pchvalues);
elseif iscell(pchvalues)
    pchvaluesCell = pchvalues;
end

% Create name-value pairs
nameValuePairs = cellfun(@(x, y) {x, y}, pchnames, pchvaluesCell, ...
                        'UniformOutput', false);

% Construct full paths to output files
[simParamsPaths, sREREsynPaths, sTCREsynPaths, sRETCsynPaths, ...
    sREspikePaths, sTCspikePaths, sREvPaths, sTCvPaths, ...
    sREcliPaths, sREsp1Paths, sREsp2Paths, sTCsp1Paths, sTCsp2Paths, ...
    sREleakPaths, sREparamsPaths, sTCparamsPaths, simCmdPaths] = ...
        argfun(@(x) cellfun(@(y) construct_fullpath(x, ...
                                    'Directory', outFolder, ...
                                    'NameValuePairs', y), ...
                            nameValuePairs, 'UniformOutput', false), ...
                simparamsF, sREREsynF, sTCREsynF, sRETCsynF, ...
                sREspikeF, sTCspikeF, sREvF, sTCvF, ...
                sREcliF, sREsp1F, sREsp2F, sTCsp1F, sTCsp2F, ...
                sREleakF, sREparamsF, sTCparamsF, scmdsF);

% If on Windows, convert all forward slashes to double forward slashes,
%   so that these paths can be used in sprintf()
if ~isunix
    [sREREsynPaths, sTCREsynPaths, sRETCsynPaths, ...
            sREspikePaths, sTCspikePaths, sREvPaths, sTCvPaths, ...
            sREcliPaths, sREsp1Paths, sREsp2Paths, sTCsp1Paths, ...
            sTCsp2Paths, sREleakPaths, sREparamsPaths, sTCparamsPaths] = ...
        argfun(@(x) replace(x, '\', '\\'), ...
                sREREsynPaths, sTCREsynPaths, sRETCsynPaths, ...
                sREspikePaths, sTCspikePaths, sREvPaths, sTCvPaths, ...
                sREcliPaths, sREsp1Paths, sREsp2Paths, sTCsp1Paths, ...
                sTCsp2Paths, sREleakPaths, sREparamsPaths, sTCparamsPaths)
end

%% Create a table for simulation parameters
simParamLabels = { ...
    '# of cells'; 'temperature of experiment (celsius)'; ...
    'number of special neurons'; 'radius of intra-RE connections'; ...
    'radius of TC-RE synaptic connections'; ...
    'radius of RE-TC synaptic connections'; ...
    'pharm condition'; 'GABA-B conductance amplitude scaling (%)'; ...
    'action potential threshold (mV)'; 'synaptic delay (ms)'; ...
    'synaptic weight (fraction of channels activated)'; ...
    'whether to use HH channels'; ...
    'number of segments in an RE cell (must be odd)'; ...
    'which cld mechanism to use (0 or 1 or 2) in RE cells'; ...
    'whether to concentrate synapses (0 or 1) in RE cells'; ...
    'Cl- removal time constant (ms) in RE cells'; ...
    'leak reversal potential (mV) of RE cells'; ...
    'diameter (um) of an RE cell'; ...
    'lower bound for passive leak conductance (S/cm^2) in RE cells'; ...
    'upper bound for passive leak conductance (S/cm^2) in RE cells'; ...
    'maximal conductance (uS) of the GABA-A receptor on RE cells'; ...
    'maximal conductance (uS) of the AMPA receptor on RE cells'; ...
    'reversal potential (mV) of the GABA-A receptor on TC cells'; ...
    'maximal conductance (uS) of the GABA-A receptor on TC cells'; ...
    'reversal potential (mV) of the GABA-B receptor on TC cells'; ...
    'conductance amplitude (uS) of the GABA-B receptor on TC cells'; ...
    'rising phase time constant (ms) of the GABA-B receptor on TC cells'; ...
    'fast decay time constant (ms) of the GABA-B receptor on TC cells'; ...
    'slow decay time constant (ms) of the GABA-B receptor on TC cells'; ...
    'weight (1) of the fast decay of the GABA-B receptor on TC cells'; ...
    'initial intracellular [Ca++] (mM)'; 'initial extracellular [Ca++] (mM)'; ...
    'initial intracellular [Cl-] (mM)'; 'initial extracellular [Cl-] (mM)'; ...
    'activation mode'; 'ID # of central neuron to activate'; ...
    'stimulation delay (ms)'; 'stimulation duration (ms)'; 'stimulation frequency (Hz)'; ...
    'current pulse duration (ms)'; 'current pulse amplitude (nA)'; ...
    'current pulse period (ms)'; 'number of current pulses'; ...
    'voltage (mV) to set activated neuron to'; ...
    'width of Gaussian distribution for randomly activating cells'; ...
    'maximum likelihood of activation at center'; ...
    'simulation mode'; 'total time of simulation (ms)'; ...
    'time step of integration (ms)'; 'time to start plotting (ms)'; ...
    'ID # of 1st special RE neuron to record'; 'ID # of 2nd special RE neuron to record'; ...
    'ID # of 1st special TC neuron to record'; 'ID # of 2nd special TC neuron to record'; ...
    'ID # of the activated neuron'; 'ID # of the neuron one below the activated neuron'; ...
    'ID # of the neuron 2 below the activated neuron'; 'ID # of a far away neuron'; ...
    'whether in debug mode'; ...
    'number to seed random number generator'; ...
    'whether to run on large memory nodes'; ...
    'whether TC neurons are heterogeneous'; ...
    'whether GABA-A conductances are removed'; ...
    'whether to save network topology'; 'whether to save spike data'; 'whether to save all voltage data'; ...
    'whether to save all chloride concentration data'; 'whether to save special neuron data'; ...
    'whether to plot network topology'; 'whether to plot spike data'; 'whether to plot single neuron data'; ...
    'number of times to run simulation'; 'current simulation number'};

% Note: Must be consistent with m3ha_network_update_dependent_params.m
simParamNames = { ...
    'nCells'; 'celsius'; 'nSpecial'; 'RERErad'; ...
    'TCRErad'; 'RETCrad'; 'pCond'; 'gIncr'; ...
    'spThr'; 'synDel'; 'synWeight'; 'useHH'; ...
    'REnsegs'; 'REcldnum'; ...
    'REconsyn'; 'REtauKCC2'; ...
    'REepas'; 'REdiam'; 'REgpasLB'; 'REgpasUB'; ...
    'REgabaGmax'; 'REampaGmax'; 'TCgabaaErev'; 'TCgabaaGmax'; 'TCgababErev'; ...
    'TCgababAmp'; 'TCgababTrise'; 'TCgababTfallFast'; 'TCgababTfallSlow'; 'TCgababW'; ...
    'cai0'; 'cao0'; 'cli0'; 'clo0'; ...
    'actMode'; 'actCellID'; ...
    'stimStart'; 'stimDur'; 'stimFreq'; ...
    'cpDur'; 'cpAmp'; 'cpPer'; 'cpNum'; ...
    'actCellV'; 'actWidth'; 'actMaxP'; ...
    'simMode'; 'tStop'; ...
    'dt'; 'tStart'; ...
    'REsp1cellID'; 'REsp2cellID'; ...
    'TCsp1cellID'; 'TCsp2cellID'; ...
    'act'; 'actLeft1'; 'actLeft2'; 'far'; ...
    'debugFlag'; 'seedNumber'; ...
    'onLargeMemFlag'; 'heteroTCFlag'; 'bicucullineFlag'; ...
    'saveNetwork'; 'saveSpikes'; 'saveSomaVoltage'; ...
    'saveSomaCli'; 'saveSpecial'; ...
    'plotNetwork'; 'plotSpikes'; 'plotSingleNeuronData'; ...
    'nSims'; 'simNumber'};

% Set initial values for all parameters
simParamsInit = [ ...
    nCells; celsius; nSpecial; RERErad; ...
    TCRErad; RETCrad; pCond; gIncr; ...
    spThr; synDel; synWeight; useHH; ...
    REnsegs; REcldnum; ...
    REconsyn; REtauKCC2; ...
    REepas; REdiam; REgpasLB; REgpasUB; ...
    REgabaGmax; REampaGmax; TCgabaaErev; TCgabaaGmax; TCgababErev; ...
    TCgababAmp; TCgababTrise; TCgababTfallFast; TCgababTfallSlow; TCgababW; ...
    cai0; cao0; cli0; clo0; ...
    actMode; actCellID; ...
    stimStart; stimDur; stimFreq; ...
    cpDur; cpAmp; cpPer; cpNum; ...
    actCellV; actWidth; actMaxP; ...
    simMode; tStop; ...
    dt; tStart; ...
    REsp1cellID; REsp2cellID; TCsp1cellID; TCsp2cellID; ...
    act; actLeft1; actLeft2; far; ...
    debugFlag; seedNumber; ...
    onLargeMemFlag; heteroTCFlag; bicucullineFlag; ...
    saveNetwork; saveSpikes; saveSomaVoltage; ...
    saveSomaCli; saveSpecial; ...
    plotNetwork; plotSpikes; plotSingleNeuronData; ...
    nSims; 0];

% Check if the labels, names and initial values have equal lengths
if numel(unique([numel(simParamLabels), ...
                numel(simParamNames), numel(simParamsInit)])) > 1
    error('simParamLabels, simParamNames and simParamsInit not equal length!');
end

% Create an initial parameters table
simTableInit = table(simParamsInit, simParamLabels, ...
                        'RowNames', simParamNames, ...
                        'VariableNames', {'Value', 'Label'});

%% Setup simulation parameters for each simulation
% Create simulation parameter tables for each simulation
simTables = ...
    cellfun(@(x, y, z) setup_params_table(simTableInit, x, y, z, ...
                                            experimentName), ...
                        pchnames, pchvaluesCell, ...
                        num2cell(transpose(1:nSims)), ...
                        'UniformOutput', false);

% Save parameter tables
cellfun(@(x, y) writetable(x, y, 'WriteRowNames', true), ...
        simTables, simParamsPaths);

%% Load cell parameters to use
% Create cell IDs for TC neurons
%   Note: Index starts from 0 in NEURON
TCcellID = transpose(1:nCells) - 1;

% Find the IDs of cells that are stimulated or artificially activated
stimCellIDs = m3ha_network_define_actmode(actMode, actCellID, nCells, ...
                                            RERErad, RETCrad);

% Seed random number generator with repetition number
rng(seedNumber);

% Generate template IDs to use for each TC neuron
templateIDsUsed = candidateIDs(randi(nCandidates, nCells, 1));

% Select templates for each TC neuron
templateNamesUsed = templateNames(templateIDsUsed);

% Create template file names
templateFileNames = strcat('bestparams_', templateNamesUsed, '.csv');

% Create full paths to the template files
templatePaths = fullfile(homeDirectory, paramsDirName, templateFileNames);

% Import NEURON parameters
TCparamTables = load_params(templatePaths);

% Create paths to NEURON parameter tables
TCparamFileNames = create_labels_from_numbers(TCcellID, 'Prefix', 'Cell', ...
                                            'Suffix', '_params.csv');
TCparamPaths = fullfile(outFolder, TCparamFileNames);

% Save NEURON parameter tables
cellfun(@(x, y) save_params(x, 'FileName', y), TCparamTables, TCparamPaths, ...
        'UniformOutput', false);

% Load parameters for TC neurons
for iCell = 1:nCells
    % Get the table for this neuron
    TCparamTableThis = TCparamTables{iCell};

    % Extract parameter names and values
    TCparamNamesThis = TCparamTableThis.Properties.RowNames;
    TCparamValuesThis = TCparamTableThis{:, 'Value'};

    % Store each parameter in an appropriately named vector
    % cellfun(@(x, y) eval(sprintf('TC%s(%d) = %g;', x, iCell, y)), ...
    %         TCparamNamesThis, num2cell(TCparamValuesThis));
    nTCParams = numel(TCparamNamesThis);
    for iParam = 1:nTCParams
        % Store parameter in the corresponding vector
        eval(sprintf('TC%s(%d) = %g;', TCparamNamesThis{iParam}, ...
                    iCell, TCparamValuesThis(iParam)));
    end
end

% Create parameters from old parameters
if ~exist('TCdiamDend', 'var') 
    TCdiamDend = TCdiamDendToSoma .* TCdiamSoma;
end

%% Build simulation commands to be read by NEURON through the here-document
% TODO: Might be too complicated to make into a function
% simCommands = cellfun(@(x, y) create_simulation_commands(x, y), ...
%                         simTables, simCmdPaths, 'UniformOutput', false);
% function simCommand = create_simulation_commands(paramsTable, scmdsPath)
simCommands = cell(nSims, 1);             % stores simulation commands
for iSim = 1:nSims
    % Retrieve simulation parameter names and values for this simulation
    simParamNames = simTables{iSim}.Properties.RowNames;
    paramValues = simTables{iSim}{:, 'Value'};

    % Update simulation parameters for this simulation
    for iParam = 1:length(paramValues)
        eval(sprintf('%s = %g;', simParamNames{iParam}, paramValues(iParam)));
    end
    simCommand = '';
    
    % Commands to create TC neurons
    for iCell = 1:nCells
        simCommand = [simCommand, sprintf(['buildTC(%d, %g, %g, %g, %g, ', ...
            '%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, ', ...
            '%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, ', ...
            '%g, %g, %g, %g, %d, %d)\n'], ...
            TCcellID(iCell), TCdiamSoma(iCell), TCLDend(iCell), ...
            TCdiamDend(iCell), TCdistDendPercent(iCell), ...
            TCcorrD(iCell), TCgpas(iCell), TCepas(iCell), ...
            TCpcabarITSoma(iCell), TCpcabarITDend1(iCell), ...
            TCpcabarITDend2(iCell), TCshiftmIT(iCell), TCshifthIT(iCell), ...
            TCslopemIT(iCell), TCslopehIT(iCell), ...
            TCghbarIhSoma(iCell), TCghbarIhDend1(iCell), ...
            TCghbarIhDend2(iCell), TCehIh(iCell), TCshiftmIh(iCell), ...
            TCgkbarIASoma(iCell), TCgkbarIADend1(iCell), ...
            TCgkbarIADend2(iCell), TCgkbarIKirSoma(iCell), ...
            TCgkbarIKirDend1(iCell), TCgkbarIKirDend2(iCell), ...
            TCgnabarINaPSoma(iCell), TCgnabarINaPDend1(iCell), ...
            TCgnabarINaPDend2(iCell), useHH, templateIDsUsed(iCell))];
    end

    % Commands to create RE neurons and build network
    simCommand = [simCommand, sprintf(['buildnet("%s", "%s", "%s", ', ...
        '%g, %g, %g, %g, %g, %d, ', ...
        '%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, ', ...
        '%g, %g, %g, %g, %g, %g, %g, %g, ', ...
        '%g, %g, %g, %g, %g, %g, %g, ', ...
        '%d, %d, %d, %d, %d, %d, %d)\n'], ...
        sREREsynPaths{iSim}, sTCREsynPaths{iSim}, sRETCsynPaths{iSim}, ...
        REsp1cellID, REsp2cellID, TCsp1cellID, TCsp2cellID, useHH, REnsegs, ...
        REcldnum, REconsyn, REtauKCC2, REepas, REdiam, ...
        REgabaGmax, REampaGmax, TCgabaaErev, TCgabaaGmax, TCgababErev, ...
        TCgababAmp, TCgababTrise, TCgababTfallFast, TCgababTfallSlow, TCgababW, ...
        RERErad, TCRErad, RETCrad, ...
        spThr, synDel, synWeight, cai0, cao0, cli0, clo0, ...
        actCellID, actMode, saveNetwork, saveSpikes, saveSomaVoltage, saveSomaCli, saveSpecial)];

    % Commands to randomize leak current properties
    %     uniformly randomizes leak conductance in [REgpasLB, REgpasUB]
    simCommand = [simCommand, sprintf('randleak(%g, %g, "%s")\n', ...
        REgpasLB, REgpasUB, sREleakPaths{iSim})];

    % Commands to initialize variables
    simCommand = [simCommand, sprintf('vinitRE(%g)\n', REepas)];

    % Commands to set up neural activation protocol
    switch actMode
    case 1
        % Activate a single RE cell by injecting a train of current pulses
        simCommand = [simCommand, sprintf('REsinglecp(%g, %g, %g, %g, %g, %g)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum)];
    case 2
        % Activate every (RERErad + 1)th RE cell by injecting trains of current pulses
        simCommand = [simCommand, sprintf('REmultcp(%g, %g, %g, %g, %g, %g, %g)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum, RERErad)];
    case 3
        % Activate 3 RE cells (RERErad + 1) apart by injecting trains of current pulses
        simCommand = [simCommand, sprintf('REthreecp(%g, %g, %g, %g, %g, %g, %g)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum, RERErad)];
    case 4
        % Activate a single RE cell at a specific voltage
        simCommand = [simCommand, sprintf('REsingleact(%g, %g)\n', ...
                    actCellID, actCellV)];    
    case 5
        % Activate RE cells with a Gaussian likelihood at a specific voltage
        simCommand = [simCommand, sprintf('RErandact(%g, %g, %g, %g)\n', ...
                    actCellID, actWidth, actMaxP, actCellV)];
    case 6
        % Activate every 3rd RE cell by injecting trains of current pulses
        simCommand = [simCommand, sprintf('REmultcp(%g, %g, %g, %g, %g, %g, %g)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum, 2)];
    case 7
        % Activate every RE cell by injecting trains of current pulses
        simCommand = [simCommand, sprintf('REmultcp(%g, %g, %g, %g, %g, %g, %g)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum, 0)];
    case 8
        % Activate 3 RE cells RETCrad apart by injecting trains of current pulses
        simCommand = [simCommand, sprintf('REthreecp(%g, %g, %g, %g, %g, %g, %g)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum, RETCrad-1)];
    case 9
        % Activate 10 center RE cells by injecting trains of current pulses
        simCommand = [simCommand, sprintf('REcentercp(%g, %g, %g, %g, %g, %g, %d)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum, 10)];
    case 10
        % Activate 20 center RE cells by injecting trains of current pulses
        simCommand = [simCommand, sprintf('REcentercp(%g, %g, %g, %g, %g, %g, %d)\n', ...
                    actCellID, stimStart, cpDur, cpAmp, cpPer, cpNum, 20)];
    otherwise
        error('actMode undefined!');
    end

    % Commands to run simulation
    %%%%%%
    %%%%%%%%%%%%
    simCommand = [simCommand, sprintf(['sim(%g, %g, "%s", "%s", "%s", "%s", "%s", ', ...
                                        '"%s", "%s", "%s", "%s", %d, %d, %d, %d)\n'], ...
                    tStop, dt, ...
                    sREspikePaths{iSim}, sTCspikePaths{iSim}, sREvPaths{iSim}, sTCvPaths{iSim}, sREcliPaths{iSim}, ...
                    sREsp1Paths{iSim}, sREsp2Paths{iSim}, sTCsp1Paths{iSim}, sTCsp2Paths{iSim}, ...
                    saveSpikes, saveSomaVoltage, saveSomaCli, saveSpecial)];
    %%%%%%%%%%%%
    %%%%%%

    % Commands to print all parameters
    simCommand = [simCommand, sprintf('print_params("%s", "%s", %g, %g, %g, %g)\n', ...
                    sREparamsPaths{iSim}, sTCparamsPaths{iSim}, useHH, REnsegs, REcldnum, REconsyn)];

    % Print simulation commands to a text file
    fid = fopen(simCmdPaths{iSim}, 'w');
    fprintf(fid, '%s\n\n', simCommand);
    fclose(fid);

    % Store in array
    simCommands{iSim} = simCommand;
end

%% Launch NEURON and execute hocFile
outputTable = run_neuron(hocFileName, 'SimCommands', simCommands, ...
                    'SimNumbers', simNumbers, 'OutFolder', outFolder, ...
                    'DebugFlag', debugFlag, 'OnHpcFlag', onHpcFlag, ...
                    'SaveStdOutFlag', saveStdOutFlag, ...
                    'RenewParpoolFlag', renewParpoolFlagNeuron, ...
                    'MaxNumWorkers', maxNumWorkersNeuron);

%% Plot stuff
timer2 = tic();

% Read data from the previous outFolder
inFolder = outFolder;

% Show network topology
if plotNetwork
    [RERE, TCRE, RETC] = m3ha_network_show_net(inFolder, 'OutFolder', outFolder, ...
                                    'FirstOnly', true);
end

% Show spike raster plot for each set of neurons (each .spi file)
[~, ~, numActive, latency, oscDur] ...
    = m3ha_network_raster_plot(inFolder, 'OutFolder', outFolder, ...
            'RenewParpool', renewParpoolFlagPlots, ...
            'MaxNumWorkers', maxNumWorkersPlots, ...
            'SingleTrialNum', simNumbers, ...
            'PlotSpikes', plotSpikes, 'PlotTuning', plotTuning);

% Show single neuron traces and heat maps for selected neurons (each .singv, .singcli & .singsp file)
if plotSingleNeuronData
    m3ha_network_single_neuron(inFolder, 'OutFolder', outFolder, ...
            'CellsToPlot', cellsToPlot, 'PropertiesToPlot', propertiesToPlot, ...
            'RenewParpool', renewParpoolFlagPlots, ...
            'MaxNumWorkers', maxNumWorkersPlots);
end

%% Compute time taken
timeTaken = toc(timer2);
fprintf('It took %3.3g seconds to plot and save stuff!!\n', timeTaken);
fprintf('\n');

%% Save all variables again in a mat file named by the date & time
if saveAllVariablesFlag
    save(fullfile(outFolder, sprintf('%s.mat', fileLabel)), '-v7.3');
end

%% Play Handel if not on Rivanna
if exist('/media/adamX/', 'dir') == 7
    load handel
    sound(y, Fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paramsTable = setup_params_table(paramsTable, pchname, pchvalue, ...
                                            simNumber, experimentName)

% Update simulation number
paramsTable{'simNumber', 'Value'} = simNumber;

% Update parameters
paramsTable = m3ha_network_change_params(paramsTable, pchname, pchvalue, ...
                                        'ExperimentName', experimentName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%% Setup parameters and filenames for each trial
paramValues = cell(1, nSims);     % stores parameters used for each trial
simParamsPaths = cell(1, nSims);
sREREsynPaths = cell(1, nSims);
sTCREsynPaths = cell(1, nSims);
sRETCsynPaths = cell(1, nSims);
sREspikePaths = cell(1, nSims);
sTCspikePaths = cell(1, nSims);
sREvPaths = cell(1, nSims);
sTCvPaths = cell(1, nSims);
sREcliPaths = cell(1, nSims);
sREsp1Paths = cell(1, nSims);
sREsp2Paths = cell(1, nSims);
sTCsp1Paths = cell(1, nSims);
sTCsp2Paths = cell(1, nSims);
sREleakPaths = cell(1, nSims);
sREparamsPaths = cell(1, nSims);
sTCparamsPaths = cell(1, nSims);
for k = 1:nSims
    % Update trial count
    idxTrialNumber = find_in_strings('simNumber', simParamNames);
    paramValues{k}(idxTrialNumber) = k;

    % Set parameters for this trial according to pchnames and pchvalues
    pchnamesThis = pchnames{k};
    if iscell(pchvalues)
        pchvaluesThis = pchvalues{k};
    elseif isnumeric(pchvalues)
        pchvaluesThis = pchvalues(k);
    end
    paramValues{k} = ...
        m3ha_network_change_params(pchnamesThis, pchvaluesThis, ...
                    simParamNames, simParamsInit, 'ExperimentName', experimentName);

    % Create name-value pairs
    nameValuePairsThis = {pchnamesThis, pchvaluesThis};

    % Print parameters to a comma-separated-value file
    simParamsPaths{k} = construct_fullpath(simparamsF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    fid = fopen(simParamsPaths{k}, 'w');
    for iParam = 1:length(paramValues{k})
        fprintf(fid, '%s, %g, %s\n', ...
            simParamNames{iParam}, paramValues{k}(iParam), simParamLabels{iParam});
    end
    fclose(fid);

    % Construct full file names
    sREREsynPaths{k}   = construct_fullpath(sREREsynF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sTCREsynPaths{k}   = construct_fullpath(sTCREsynF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sRETCsynPaths{k}   = construct_fullpath(sRETCsynF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sREspikePaths{k}   = construct_fullpath(sREspikeF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sTCspikePaths{k}   = construct_fullpath(sTCspikeF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sREvPaths{k}       = construct_fullpath(sREvF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sTCvPaths{k}       = construct_fullpath(sTCvF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sREcliPaths{k}     = construct_fullpath(sREcliF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sREsp1Paths{k}     = construct_fullpath(sREsp1F, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sREsp2Paths{k}     = construct_fullpath(sREsp2F, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sTCsp1Paths{k}     = construct_fullpath(sTCsp1F, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sTCsp2Paths{k}     = construct_fullpath(sTCsp2F, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sREleakPaths{k}     = construct_fullpath(sREleakF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sREparamsPaths{k}     = construct_fullpath(sREparamsF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
    sTCparamsPaths{k}     = construct_fullpath(sTCparamsF, 'Directory', outFolder, ...
                        'NameValuePairs', nameValuePairsThis);
end

% Import .p file data
pfiledata = importdata(sprintf('pfiles/bestparams_%s.p', ...
                                templateNamesUsed{iCell}));  % a XX x 4 array
nTCParams = size(pfiledata, 1);
for k = 1:nTCParams                 % for each parameter saved
    % Store parameter in the corresponding vector
    eval(sprintf('TC%s(%d) = %g;', pfiledata{k, 1}, iCell, pfiledata{k, 2}));
end

paramsIn = paramsTable{:, 'Value'};
simParamNames = paramsTable.Properties.RowNames;
paramsOut = m3ha_network_change_params(pchname, pchvalue, ...
                    simParamNames, paramsIn, 'ExperimentName', experimentName);
paramsTable{:, 'Value'} = paramsOut;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
