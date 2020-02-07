function [REspikes, TCspikes, numActiveTC, numActiveRE, ...
            latencyTC, latencyRE, oscDurTC, oscDurRE, nSpikesTC, nSpikesRE, ...
            actDurTC, actDurRE, actVelTC, actVelRE] = ...
                m3ha_network_raster_plot(inFolder, varargin)
%% Shows a spike raster plot and compute numActive, latency, oscDur, nSpikes, actDur, actVel for each set of neurons (each .spi file in the inFolder)
% USAGE: [REspikes, TCspikes, numActiveTC, numActiveRE, ...
%           latencyTC, latencyRE, oscDurTC, oscDurRE, nSpikesTC, nSpikesRE, ...
%           actDurTC, actDurRE, actVelTC, actVelRE] = ...
%               m3ha_network_raster_plot(inFolder, varargin)
% Arguments:
%    inFolder   - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%               must be a directory
%    varargin   - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%               could be anything recognised by the built-in saveas() function
%               (see isfigtype.m under Adams_Functions)
%               default == 'png'
%               - 'OutFolder': the name of the directory that the plots will be placed
%               must be a directory
%               default: same as inFolder
%               - 'MaxNumWorkers': maximum number of workers for running NEURON 
%               must a positive integer
%               default: 20
%               - 'SingleTrialNum': number of single trial ran
%                                       or 0 if not just a single trial
%               must a positive integer
%               default: 0
%               - 'RenewParpool': whether to renew parpool every batch to release memory
%               must be logical 1 (true) or 0 (false)
%               default: true, but alway false if plotSpikes == false
%               - 'PlotSpikes': whether to plot raster plots
%               must be logical 1 (true) or 0 (false)
%               default == true
%               - 'PlotTuning': whether to plot tuning curves
%               must be logical 1 (true) or 0 (false)
%               default == true
%               - 'PlotOnly': whether to plot only and not create the figure
%               must be logical 1 (true) or 0 (false)
%               default == false
%               - 'RastersOnly': whether to plot rasters only
%               must be logical 1 (true) or 0 (false)
%               default == false
%               - 'NoRasters': whether to omit plotting rasters
%               must be logical 1 (true) or 0 (false)
%               default == false
%
% Requires:
%       cd/argfun.m
%       cd/combine_strings.m
%       cd/create_subplots.m
%       cd/find_in_strings.m
%       cd/extract_looped_params.m
%       cd/isfigtype.m
%       cd/m3ha_network_define_actmode.m
%       cd/m3ha_network_tuning_curves.m
%       cd/m3ha_network_tuning_maps.m
%       cd/plot_window_boundaries.m
%       cd/save_all_figtypes.m
%       inFolder/*.spi
%       inFolder/*loopedparams.mat
%       inFolder/['sim_params_', pstring, '.csv'] for all the possible parameter strings
%
% Used by:
%       cd/m3ha_network_launch.m
%       cd/m3ha_network_tuning_curves.m
%       cd/m3ha_network_tuning_maps.m

% File History:
% 2017-10-31 Modified from /RTCl/raster_plot.m
% 2017-11-03 Removed activation velocity and added oscillation oscDur
% 2017-11-06 Now uses define_actmode.m
% 2017-11-07 Fixed oscillation duration
% 2018-01-24 BT - Added nSpikes and actDur
% 2018-01-26 BT - Differentiate between TC and RE cells
% 2018-01-31 BT - Added actVel, TC/RE for numActive and oscDur
% 2020-02-06 Added 'PlotOnly', 'RastersOnly', 'NoRasters' as optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'inFolder', @isdir);            % the name of the directory containing the .syn files

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigTypes', 'png', ...         % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'OutFolder', '@inFolder', @isdir); % the name of the directory that the plots will be placed
addParameter(iP, 'MaxNumWorkers', 20, ...       % maximum number of workers for running NEURON
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'RenewParpool', false, ...      % whether to renew parpool every batch to release memory
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SingleTrialNum', [], ...      % number of single trial ran
    @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));
addParameter(iP, 'PlotSpikes', true, ...        % whether to plot raster plots
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotTuning', true, ...        % whether to plot tuning curves
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotOnly', false, ...         % whether to just plot
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RastersOnly', false, ...         % whether to just plot
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'NoRasters', false, ...         % whether to just plot
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, inFolder, varargin{:});
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
outFolder = iP.Results.OutFolder;
maxNumWorkers = iP.Results.MaxNumWorkers;
renewParpool = iP.Results.RenewParpool;
singleTrialNum = iP.Results.SingleTrialNum;
plotSpikes = iP.Results.PlotSpikes;
plotTuning = iP.Results.PlotTuning;
plotOnly = iP.Results.PlotOnly;
rastersOnly = iP.Results.RastersOnly;
noRasters = iP.Results.NoRasters;

% Change default arguments if necessary
if strcmp(outFolder, '@inFolder')
    outFolder = inFolder;
end
if ~plotSpikes || (numel(singleTrialNum) >= 1 && singleTrialNum ~= 0)
    renewParpool = false;        % only plotting spikes will take up significant amount of memory
end

%% Find all .spi files
REfiles = dir(fullfile(inFolder, 'RE*.spi'));
nREfiles = length(REfiles);
TCfiles = dir(fullfile(inFolder, 'TC*.spi'));
nTCfiles = length(TCfiles);

%% Get loop parameters
[nump, pnames, plabels, pislog, pvalues, nperp, pchnames, pchvalues, nCells, actMode, loopmode] = extract_looped_params(inFolder);

%% Check number of .spi files
nSims = numel(pchnames);
if nREfiles ~= nSims
    fprintf('Warning: Number of RE*.spi files and number of trials don''t match!\n');
end
if nTCfiles ~= nSims
    fprintf('Warning: Number of TC*.spi files and number of trials don''t match!\n');
end

%% Create raster plot
REspikes = cell(nSims, 1);        % raw spike data for RE neurons
TCspikes = cell(nSims, 1);        % raw spike data for TC neurons
latencyTC = zeros(nSims, nCells); % latency to activation (seconds) for each TC cell
latencyRE = zeros(nSims, nCells); % latency to activation (seconds) for each RE cell
nSpikesTC = zeros(nSims, nCells); % number of spikes for each TC cell
nSpikesRE = zeros(nSims, nCells); % number of spikes for each RE cell
actDurTC = zeros(nSims, nCells);  % active duration (seconds) for each TC cell
actDurRE = zeros(nSims, nCells);  % active duration (seconds) for each RE cell
numActiveTC = zeros(nSims, 1);    % number of cells activated
numActiveRE = zeros(nSims, 1);	% number of cells activated
oscDurTC = zeros(nSims, 1);       % oscillation duration (seconds)
oscDurRE = zeros(nSims, 1);       % oscillation duration (seconds)
actVelTC = zeros(nSims, 1);       % activation velocity (cells/seconds)
actVelRE = zeros(nSims, 1);       % activation velocity (cells/seconds)
ct = 0;                             % counts number of trials completed
poolObj = gcp('nocreate');          % get current parallel pool object without creating a new one
if isempty(poolObj)
    poolObj = parpool;              % create a default parallel pool object
    oldNumWorkers = poolObj.NumWorkers;    % number of workers in the default parallel pool object
else
    oldNumWorkers = poolObj.NumWorkers;    % number of workers in the current parallel pool object
end
numWorkers = min(oldNumWorkers, maxNumWorkers);    % number of workers to use for running NEURON
if renewParpool
    delete(poolObj);                % delete the parallel pool object to release memory
end
while ct < nSims                  % while not trials are completed yet
    if ~isempty(singleTrialNum) && singleTrialNum ~= 0
        first = singleTrialNum; % run that trial
    else
        first = ct + 1;         % first trial in this batch
    end
    if ~isempty(singleTrialNum) && singleTrialNum ~= 0
        last = singleTrialNum;      % run only that trial
    elseif renewParpool && ...
        ct + numWorkers <= nSims  % if memory is to be released
        last = ct + numWorkers;     % limit the batch to numWorkers
    else
        last = nSims;             % run all trials at once
    end
    if renewParpool
        % Recreate a parallel pool object using fewer workers to prevent running out of memory
        poolObj = parpool('local', numWorkers);    
    end
    %parfor i = first:last
    for i = first:last
        % Construct current parameter string
        pstring = '';               % initialize for parfor
        if iscell(pchvalues)
            pstring = combine_strings('NameValuePairs', {pchnames{i}, pchvalues{i}});
        elseif isnumeric(pchvalues)
            pstring = combine_strings('NameValuePairs', {pchnames{i}, pchvalues(i)});
        end
        pspifile = [pstring, '.spi'];

        % Construct full path to simulation parameters spreadsheet
        simParamsPath = fullfile(inFolder, ['sim_params_', pstring, '.csv']);

        % Load spike data
        jnowRE = 0;                 % initialize for parfor
        jnowTC = 0;                 % initialize for parfor
        for j = 1:nREfiles
            if ~isempty(strfind(REfiles(j).name, pspifile))
                jnowRE = j;         % index in REfiles of current file
                REspikes{i} = load(fullfile(inFolder, REfiles(jnowRE).name));
            end
        end
        for j = 1:nTCfiles
            if ~isempty(strfind(TCfiles(j).name, pspifile))
                jnowTC = j;         % index in TCfiles of current file
                TCspikes{i} = load(fullfile(inFolder, TCfiles(jnowTC).name));
            end
        end

        % Don't plot if file not found or no spikes
        stimcellIDs = 0;            % initialize for parfor    
        if isempty(REspikes{i}) && isempty(TCspikes{i})
            fprintf('The file %s does not exist or has no spikes!\n', pspifile);
            numActiveTC(i) = 0;     % number of cells activated is zero
			numActiveRE(i) = 0;		% number of cells activated is zero
            latencyTC(i, :) = Inf(1, nCells);    
                                    % latency to activation (seconds) is infinite for every cell
            latencyRE(i, :) = Inf(1, nCells);    
                                    % latency to activation (seconds) is infinite for every cell
            oscDurTC(i) = 0;        % oscillation duration (seconds) is zero
            oscDurRE(i) = 0;        % oscillation duration (seconds) is zero
			nSpikesTC(i, :) = Inf(1, nCells); % number of spikes for each cell
			nSpikesRE(i, :) = Inf(1, nCells); % number of spikes for each cell
			actDurTC(i, :) = Inf(1, nCells);	% active duration (seconds) for each cell
			actDurRE(i, :) = Inf(1, nCells);	% active duration (seconds) for each cell
			actVelTC(i) = 0; 	% activation velocity (cells/seconds)
			actVelRE(i) = 0; 	% activation velocity (cells/seconds)
        else
            % Read the simulation parameters table
            simParamsTable = readtable(simParamsPath, 'ReadRowNames', true);

            % Extract parameters from the simulation parameters table
            [tStart, tStop, RERErad, RETCrad, ...
                actCellID, stimStart, stimDur, stimFreq] = ...
                    argfun(@(x) simParamsTable{x, 'Value'}, ...
                            'tStart', 'tStop', 'RERErad', 'RETCrad', ...
                            'actCellID', 'stimStart', 'stimDur', 'stimFreq');

            % Extract info
            if ~isempty(REspikes{i})
                REspikecelln = REspikes{i}(:, 1);
                REspiketimes = REspikes{i}(:, 2);
            else
                REspikecelln = [];
                REspiketimes = [];
            end
            if ~isempty(TCspikes{i})
                TCspikecelln = TCspikes{i}(:, 1);
                TCspiketimes = TCspikes{i}(:, 2);
            else
                TCspikecelln = [];
                TCspiketimes = [];
            end
            xLimitsLeft = tStart;
            xLimitsRight = tStop;
            stimStartPlot = stimStart;
            stimDurPlot = stimDur;
            xLabel = 'Spike time (ms)';
    
            % Find number of cells activated
            %   i.e., number of unique cells with spike times
			ind = find(TCspiketimes > stimStart);
			numActiveTC(i) = length(unique(TCspikecelln(ind)));
			ind = find(REspiketimes > stimStart);
			numActiveRE(i) = length(unique(REspikecelln(ind)));  

            % Find latency to activation (seconds) for each TC and RE cell
            latencyTCThis = zeros(1, nCells);    % for parfor
            for j = 1:nCells        % for each cellID == j-1
                ind = find(TCspikecelln == j-1);
                if isempty(ind)    
                    % Make latency infinite if cell not activated
                    latencyTCThis(j) = Inf;
                else
                    latencyTCThis(j) = (min(TCspiketimes(ind)) - stimStart) / 1000;
                end
            end
            latencyTC(i, :) = latencyTCThis;       % for parfor

            latencyREThis = zeros(1, nCells);    % for parfor
            for j = 1:nCells        % for each cellID == j-1
                ind = find(REspikecelln == j-1);
                if isempty(ind)    
                    % Make latency infinite if cell not activated
                    latencyREThis(j) = Inf;
                else
                    latencyREThis(j) = (min(REspiketimes(ind)) - stimStart) / 1000;
                end
            end
            latencyRE(i, :) = latencyREThis;       % for parfor

            % Find oscillation duration for each TC and RE cell (seconds)
            %   i.e., maximum time that any neuron has spiked minus the start of stimulation
            if isempty(TCspiketimes) || max(TCspiketimes) - stimStart < 0
                oscDurTC(i) = 0;
            else
                oscDurTC(i) = (max(TCspiketimes) - stimStart) / 1000;
			end
			if isempty(REspiketimes) || max(REspiketimes) - stimStart < 0
			    oscDurRE(i) = 0;
		    else
			    oscDurRE(i) = (max(REspiketimes) - stimStart) / 1000;
		    end

			% Find number of spikes for each TC and RE cell
            nSpikesTCThis = zeros(1, nCells);    % for parfor
            for j = 1:nCells        % for each cellID == j-1
                ind = find(TCspikecelln == j-1);
                if isempty(ind)    
                    % Make nSpikes 0 if cell not activated
                    nSpikesTCThis(j) = 0;
                else
					nSpikesTCThis(j) = length(ind);
                end
            end
            nSpikesTC(i, :) = nSpikesTCThis;       % for parfor		

            nSpikesREThis = zeros(1, nCells);    % for parfor
            for j = 1:nCells        % for each cellID == j-1
                ind = find(REspikecelln == j-1);
                if isempty(ind)    
                    % Make nSpikes 0 if cell not activated
                    nSpikesREThis(j) = 0;
                else
					nSpikesREThis(j) = length(ind);
                end
            end
            nSpikesRE(i, :) = nSpikesREThis;       % for parfor	

			% Find active duration (seconds) for each TC and RE cell
            actDurTCThis = zeros(1, nCells);    % for parfor
            for j = 1:nCells        % for each cellID == j-1
                ind = find(TCspikecelln == j-1);
                if isempty(ind)    
                    % Make actDur 0 if cell not activated
                    actDurTCThis(j) = 0;
                else
					actDurTCThis(j) = (max(TCspiketimes(ind)) - stimStart) / 1000;
                end
            end
            actDurTC(i, :) = actDurTCThis;       % for parfor	

            actDurREThis = zeros(1, nCells);    % for parfor
            for j = 1:nCells        % for each cellID == j-1
                ind = find(REspikecelln == j-1);
                if isempty(ind)    
                    % Make actDur 0 if cell not activated
                    actDurREThis(j) = 0;
                else
					actDurREThis(j) = (max(REspiketimes(ind)) - stimStart) / 1000;
                end
            end
            actDurRE(i, :) = actDurREThis;       % for parfor	

			% Find activation velocity (cells/seconds)
			if numActiveTC(i) == 0
			    actVelTC(i) = 0;
	        else
	            indStimDurCellTimes = find((TCspiketimes < stimStart + stimDur)&(TCspiketimes > stimStart));
	            stim_dur_celln = TCspikecelln(indStimDurCellTimes);
	            dn = numActiveTC(i) - length(unique(stim_dur_celln));	% difference between active cells and those stimulated at start
	            minspiketime = Inf;
	            for celln = unique(TCspikecelln)'
	                cellnind = find(TCspikecelln == celln);
	                celltimes = TCspiketimes(cellnind);
	                celltimes(celltimes <= stimStart + stimDur) = [];
	                if min(celltimes) < minspiketime
	                    minspiketime = min(celltimes);
                    end
                end
	            dt = (minspiketime - stimStart + stimDur) / 1000;
                actVelTC(i) = dn / dt;
			end
			if numActiveRE(i) == 0
			    actVelRE(i) = 0;
		    else
	            indStimDurCellTimes = find((REspiketimes < stimStart + stimDur)&(REspiketimes > stimStart));
	            stim_dur_celln = REspikecelln(indStimDurCellTimes);
	            dn = numActiveRE(i) - length(unique(stim_dur_celln));	% difference between active cells and those stimulated at start
	            minspiketime = Inf;
	            for celln = unique(REspikecelln)'
	                cellnind = find(REspikecelln == celln);
	                celltimes = REspiketimes(cellnind);
	                celltimes(celltimes <= stimStart + stimDur) = [];
	                if min(celltimes) < minspiketime
	                    minspiketime = min(celltimes);
                    end
                end
	            dt = (minspiketime - stimStart + stimDur) / 1000;
                actVelRE(i) = dn / dt;
            end
            
            % Change units of time axis from ms to s if the total time is > 10 seconds
            if tStop > 10000
                REspiketimes = REspiketimes/1000;
                TCspiketimes = TCspiketimes/1000;
                xLabel = 'Spike time (s)';
                xLimitsLeft = xLimitsLeft/1000;
                xLimitsRight = xLimitsRight/1000;
                stimStartPlot = stimStart/1000;
                stimDurPlot = stimDur/1000;
            end

            % Find maximum and minimum time to plot
            minSpikeTime = min([REspiketimes; TCspiketimes]);
            maxSpikeTime = max([REspiketimes; TCspiketimes]);
            xLimits = [min([xLimitsLeft, minSpikeTime]), ...
                        max([xLimitsRight, maxSpikeTime])];

            % Find the IDs of cells that are stimulated or artificially activated
            stimcellIDs = m3ha_network_define_actmode (actMode, actCellID, nCells, RERErad, RETCrad);

            % Separate the spike times and cell numbers into two vectors
            Lia1 = ismember(REspikecelln, stimcellIDs);    % whether the index in spikecelln belongs to stimcellIDs
            Lia2 = (REspiketimes <= stimStartPlot + stimDurPlot) ...
                & (REspiketimes >= stimStartPlot);        % whether the spike time is within the stimulation period
            REspiketimesStim = REspiketimes(Lia1 & Lia2);
            REspiketimesNonstim = REspiketimes(~(Lia1 & Lia2));
            REspikecellnStim = REspikecelln(Lia1 & Lia2);
            REspikecellnNonstim = REspikecelln(~(Lia1 & Lia2));

            % Plot raster plot
            if plotSpikes
                % Decide on the cell the plot
                if actCellID > 0
                    cellIDToPlot = actCellID - 1;
                else
                    cellIDToPlot = actCellID;
                end

                % Find the stimulation window
                stimWindow = [stimStartPlot, stimStartPlot + stimDurPlot];

                % Create plot
                if ~plotOnly
                    figNumber = 30000;
                    clearFigure = true;
                else
                    figNumber = [];
                    clearFigure = false;
                end

                % Create subplots without expanding
                [fig, ax] = create_subplots(2, 1, 'FigNumber', figNumber, ...
                                            'ClearFigure', clearFigure, ...
                                            'FigExpansion', [1, 1], ...
                                            'ExpandFromDefault', false);

                % Plot RT spikes
                plot_spikes_one_layer(ax(1), 'RT Layer', xLimits, xLabel, ...
                                REspikecellnNonstim, REspiketimesNonstim, ...
                                REspikecellnStim, REspiketimesStim, ...
                                nCells, stimWindow, stimFreq, ...
                                rastersOnly, noRasters);

                % Plot TC spikes
                plot_spikes_one_layer(ax(2), 'TC Layer', xLimits, xLabel, ...
                                TCspikecelln, TCspiketimes, ...
                                [], [], ...
                                nCells, stimWindow, [], ...
                                rastersOnly, noRasters);
%{
                % RE text labels
                text(xLimitsLeft + 0.5, nCells*1.925, ['Number of cells activated: ', num2str(numActiveRE(i))], ...
                    'Color', 'g');                          % text for number of RE cells activated
                text(xLimitsLeft + 0.5, nCells*1.875, ...
                    ['Latency to activation for Cell#', num2str(cellIDToPlot), ': ', ...
                    num2str(latencyREThis(cellIDToPlot+1)), ' seconds'], ...
                    'Color', 'g');                          % text for RE latency to activation (seconds)
                text(xLimitsLeft + 0.5, nCells*1.825, ...
                    ['Oscillation Duration: ', num2str(oscDurRE(i)), ' seconds'], ...
                    'Color', 'g');                          % text for RE oscillation duration (seconds)
                text(xLimitsLeft + 0.5, nCells*1.775, ...
                    ['Number of spikes for Cell#', num2str(cellIDToPlot), ': ', ...
                    num2str(nSpikesREThis(cellIDToPlot+1))], ...
                    'Color', 'g');                          % text for RE number of spikes 
                text(xLimitsLeft + 0.5, nCells*1.725, ...
                    ['Activation duration for Cell#', num2str(cellIDToPlot), ': ', ...
                    num2str(actDurREThis(cellIDToPlot+1)), ' seconds'], ...
                    'Color', 'g');                          % text for RE activation duration (seconds)
                text(xLimitsLeft + 0.5, nCells*1.675, ...
                    ['Activation velocity: ', num2str(actVelRE(i)), ' cells/second'], ...
                    'Color', 'g');                          % text for RE activation velocity (cells/second) 

                % TC text labels
                text(xLimitsLeft + 0.5, nCells*0.9625, ['Number of cells activated: ', num2str(numActiveTC(i))], ...
                    'Color', 'g');                          % text for number of TC cells activated
                text(xLimitsLeft + 0.5, nCells*0.9125, ...
                    ['Latency to activation for Cell#', num2str(cellIDToPlot), ': ', ...
                    num2str(latencyTCThis(cellIDToPlot+1)), ' seconds'], ...
                    'Color', 'g');                          % text for TC latency to activation (seconds)
                text(xLimitsLeft + 0.5, nCells*0.8625, ...
                    ['Oscillation Duration: ', num2str(oscDurTC(i)), ' seconds'], ...
                    'Color', 'g');                          % text for TC oscillation duration (seconds)
                text(xLimitsLeft + 0.5, nCells*0.8125, ...
                    ['Number of spikes for Cell#', num2str(cellIDToPlot), ': ', ...
                    num2str(nSpikesTCThis(cellIDToPlot+1))], ...
                    'Color', 'g');                          % text for TC number of spikes 
                text(xLimitsLeft + 0.5, nCells*0.7625, ...
                    ['Activation duration for Cell#', num2str(cellIDToPlot), ': ', ...
                    num2str(actDurTCThis(cellIDToPlot+1)), ' seconds'], ...
                    'Color', 'g');                          % text for TC activation duration (seconds)
                text(xLimitsLeft + 0.5, nCells*0.7125, ...
                    ['Activation velocity: ', num2str(actVelTC(i)), ' cells/second'], ...
                    'Color', 'g');                          % text for TC activation velocity (cells/second) 
%}

                if ~plotOnly && ~rastersOnly
                    figbaseWext = replace(REfiles(jnowRE).name, 'RE_', '');
                    figbase = replace(figbaseWext, '.spi', '');
                    suptitle(replace(figbase, '_', '\_'));
                end

                if ~plotOnly
                    % Save figure
                    figname = fullfile(outFolder, replace(figbaseWext, '.spi', '_raster_plot.png'));
                    save_all_figtypes(fig, figname, figtypes);
                    % close(fig);
                end
            end    
        end
%        close all;
    end

    if renewParpool
        delete(poolObj);    % delete the parallel pool object to release memory
    end
    if ~isempty(singleTrialNum) && singleTrialNum ~= 0
        ct = nSims + 1;   % don't run again
    else
        ct = last;          % update number of trials completed
    end
end
if renewParpool
    % Recreate a parallel pool object using the previous number of workers
    poolObj = parpool('local', oldNumWorkers);
end

%% Plot tuning curves
if plotTuning
    switch loopmode
    case 'cross'
        cellsToPlot = 0:nCells/2-1;        % cells to plot for the latency tuning curve
        m3ha_network_tuning_curves (inFolder, outFolder, numActiveTC, numActiveRE, ...
				latencyTC, latencyRE, oscDurTC, oscDurRE, nSpikesTC, nSpikesRE, ...
				actDurTC, actDurRE, actVelTC, actVelRE, ...
                nump, pnames, plabels, pislog, pvalues, nperp, ...
                nCells, actMode, loopmode, 'CellsToPlot', cellsToPlot, ...
                'FigTypes', figtypes);
    case 'grid'
        m3ha_network_tuning_maps (inFolder, outFolder, numActiveTC, numActiveRE, ...
				latencyTC, latencyRE, oscDurTC, oscDurRE, nSpikesTC, nSpikesRE, ...
				actDurTC, actDurRE, actVelTC, actVelRE, ...
                nump, pnames, plabels, pislog, pvalues, nperp, ...
                nCells, actMode, loopmode, 'FigTypes', figtypes);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_spikes_one_layer(axHandle, subTitle, xLimits, xLabel, ...
                                cellNumNonStim, spikeTimesNonstim, ...
                                cellNumStim, spikeTimesStim, ...
                                nCells, stimWindow, stimFreq, ...
                                rastersOnly, noRasters)

% Hard-coded parameters
% Pull up subplot
subplot(axHandle); hold on;

% Compute appropriate y axis limits
yLimits = [-0.5, nCells - 0.5];
yTicks = [0, round(nCells/4), round(nCells/2), round(nCells*3/4), nCells - 1];

% Plot stimulation window
if ~isempty(stimWindow) && ~noRasters
    plot_window_boundaries(stimWindow, 'YLimits', yLimits, ...
                            'BoundaryType', 'verticalLines', ...
                            'LineStyle', '-', 'Color', 'PaleGreen');
end

% Plot stimulation info text
if ~isempty(stimFreq) && ~rastersOnly
    text(stimWindow(2) + 0.5, nCells * 0.975, ...
        ['Stim: ', num2str(stimFreq), ' Hz'], 'Color', 'r');
end

% Plot spikes from non-stimulated cells
if ~isempty(spikeTimesNonstim) && ~noRasters
    line([spikeTimesNonstim, spikeTimesNonstim]', ...
        [cellNumNonStim - 0.5, cellNumNonStim + 0.5]', ...
        'Color', 'k', 'LineWidth', 0.5);
end

% Plot spikes from stimulated cells
if ~isempty(spikeTimesStim) && ~noRasters
    line([spikeTimesStim, spikeTimesStim]', ...
        [cellNumStim - 0.5, cellNumStim + 0.5]', ...
        'Color', 'r', 'LineWidth', 0.5);
end

% Set x axis limits
xlim(xLimits);

% Set y axis limits
ylim(yLimits);

% Set y tick locations
yticks(yTicks);

% Plot axis labels
if ~rastersOnly
    xlabel(xLabel);
    ylabel('Neuron ID');
end

% Show a title for the subplot
if ~rastersOnly
    title(subTitle);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
