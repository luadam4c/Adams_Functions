function [pvalues, numActiveTC, numActiveRE, latencyTC, latencyRE, ...
            oscDurTC, oscDurRE, nSpikesTC, nSpikesRE, actDurTC, actDurRE, ...
            actVelTC, actVelRE, nump, pnames, plabels, pislog, ...
            nperp, ncells, actmode, latency_cells_to_plot] = ...
                m3ha_network_tuning_curves (infolders, outfolder, numActiveTC, ...
                    numActiveRE, latencyTC, latencyRE, oscDurTC, oscDurRE, ...
                    nSpikesTC, nSpikesRE, actDurTC, actDurRE, actVelTC, ...
                    actVelRE, nump, pnames, plabels, pislog, pvalues, nperp, ...
                    ncells, actmode, loopmode, varargin)
%% Shows a tuning curve for numActive, latency, oscDur, for each parameter changed
% USAGE: [pvalues, numActiveTC, numActiveRE, latencyTC, latencyRE, ...
%           oscDurTC, oscDurRE, nSpikesTC, nSpikesRE, actDurTC, actDurRE, ...
%           actVelTC, actVelRE, nump, pnames, plabels, pislog, ...
%           nperp, ncells, actmode, latency_cells_to_plot] = ...
%               m3ha_network_tuning_curves (infolders, outfolder, numActiveTC, ...
%                   numActiveRE, latencyTC, latencyRE, oscDurTC, oscDurRE, ...
%                   nSpikesTC, nSpikesRE, actDurTC, actDurRE, actVelTC, ...
%                   actVelRE, nump, pnames, plabels, pislog, pvalues, nperp, ...
%                   ncells, actmode, loopmode, varargin)
% Arguments:
%       infolders   - the directory or a cell array of directories that contains the data to plot
%       outfolder   - (opt) the directory to place the output plots
%       TODO: other optional arguments; TODO: argument checks
%       varargin    - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - 'CellsToPlot': cell(s) to plot for latencies
%                   must be a nonnegative integer vector
%                   default == 0:(ncells/2-1) if infolder is a cell array
%                   otherwise == ncells/2
%                   - 'PCondNames': names of parameters to be varied across conditions, e.g., {'stim_freq', 'tau'}
%                   must be a cell array of strings or character arrays
%                   default == {}
%                   NOTE: If infolders is a cell array, 
%                   'PCondNames' must be provided unless the same two parameters are varied inside each infolder
%
% Requires:
%        cd/check_and_collapse_identical_contents.m
%        cd/find_in_strings.m
%        cd/get_loopedparams.m
%        cd/isfigtype.m
%        cd/m3ha_network_raster_plot.m
%        cd/plot_tuning_curve.m
%        cd/save_all_figtypes.m
%        cd/vec2cell.m
%        infolder/*loopedparams.mat
%        /home/Matlab/Downloaded_Functions/dirr.m
%
% Used by:
%        cd/m3ha_network_raster_plot.m

% File History:
% 2017-10-23 Modified from /RTCl/tuning_curves.m
% 2017-11-01 ylimits for numActive now goes to 2*ncells
% 2018-01-31 BT - Updated with nSpikes, actDur, actVel and TC/RE for all
%

%% Set readout labels
numActiveTC_label = 'Number of activated TC cells';
latencyTC_label = 'TC Latency to activation (seconds)';
oscDurTC_label = 'TC Oscillation duration (seconds)';
nSpikesTC_label = 'TC Number of spikes';
actDurTC_label = 'TC Active duration (seconds)';
actVelTC_label = 'TC Activation velocity (cells/second)';
numActiveRE_label = 'Number of activated RE cells';
latencyRE_label = 'RE Latency to activation (seconds)';
oscDurRE_label = 'RE Oscillation duration (seconds)';
nSpikesRE_label = 'RE Number of spikes';
actDurRE_label = 'RE Active duration (seconds)';
actVelRE_label = 'RE Activation velocity (cells/second)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));    % for dirr.m
addpath(fullfile(functionsdirectory, '/Adams_Functions/'));        % for isfigtype.m, save_all_figtypes.m,
                                    %    find_in_strings.m, vec2cell.m
                                    %     check_and_collapse_identical_contents.m 
                                    %    & plot_tuning_curve.m & get_loopedparams.m

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error('Not enough input arguments, type ''help m3ha_network_tuning_curves'' for usage');
end

%% Check to see if infolders is a cell array; if not, make it a cell array
if iscell(infolders)
    across_conditions = 1;      % whether the tuning curves across conditions are compared
    ncond = numel(infolders);   % number of conditions to compare
else
    infolders = {infolders};
    across_conditions = 0;
end
%% Set outfolder and create if doesn't exist
if nargin < 2 || isempty(outfolder)
    if across_conditions
        outfolder = 'tuning_curve_across_conditions';
    else
        outfolder = infolders{1};
    end
end
if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
end
%% Get numActive, latency, oscDur if not provided
if isempty(oscDurTC)
    if across_conditions
        numActiveTC = cell(1, ncond);
        numActiveRE = cell(1, ncond);
        latencyTC = cell(1, ncond);
        latencyRE = cell(1, ncond);
        oscDurTC = cell(1, ncond);
        oscDurRE = cell(1, ncond);
        nSpikesTC = cell(1, ncond);
        nSpikesRE = cell(1, ncond);
        actDurTC = cell(1, ncond);
        actDurRE = cell(1, ncond);
        actVelTC = cell(1, ncond);
        actVelRE = cell(1, ncond);
        for k = 1:ncond
            [~, ~, numActiveTC{k}, numActiveRE{k}, latencyTC{k}, latencyRE{k}, oscDurTC{k}, oscDurRE{k}, ...
                nSpikesTC{k}, nSpikesRE{k}, actDurTC{k}, actDurRE{k}, actVelTC{k}, actVelRE{k}] = ...
                    raster_plot(infolders{k}, 'OutFolder', outfolder, ...
                        'RenewParpool', 0, 'PlotSpikes', 0, 'PlotTuning', 0);
        end
    else
        [~, ~, numActiveTC, numActiveRE, latencyTC, latencyRE, oscDurTC, oscDurRE, ...
            nSpikesTC, nSpikesRE, actDurTC, actDurRE, actVelTC, actVelRE] = ...
                raster_plot(infolders{1}, 'OutFolder', outfolder, ...
                    'RenewParpool', 0, 'PlotSpikes', 0, 'PlotTuning', 0);
    end
end

%% Get looped parameters info if not provided
if isempty(actmode)        %%% TODO: might need more isempty()
    if across_conditions
        % Invariants for each condition
        nump = zeros(1, ncond);
        pnames = cell(1, ncond);
        plabels = cell(1, ncond);
        pislog = cell(1, ncond);
        ncells = zeros(1, ncond);
        actmode = zeros(1, ncond);    
        loopmode = cell(1, ncond);    

        % Variants for each condition
        pvalues = cell(1, ncond);   % parameter values might be different for each condition
        nperp = cell(1, ncond);     % number of parameter values might be different for each condition

        % Get the parameters from each condition
        for k = 1:ncond
            [nump(k), pnames{k}, plabels{k}, pislog{k}, pvalues{k}, nperp{k}, ...
                ~, ~, ncells(k), actmode(k), loopmode{k}] = get_loopedparams(infolders{k});
        end
    else
        [nump, pnames, plabels, pislog, pvalues, nperp, ~, ~, ncells, actmode, loopmode] = get_loopedparams(infolders{1});
    end
end
% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'infolders', ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));            % the name of the directory containing the .syn files
addRequired(iP, 'outfolder', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(iP, 'numActiveTC', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'numActiveRE', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'latencyTC', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', '2d'}));
addRequired(iP, 'latencyRE', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', '2d'}));
addRequired(iP, 'oscDurTC', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'oscDurRE', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'nSpikesTC', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', '2d'}));
addRequired(iP, 'nSpikesRE', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', '2d'}));
addRequired(iP, 'actDurTC', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', '2d'}));
addRequired(iP, 'actDurRE', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', '2d'}));
addRequired(iP, 'actVelTC', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'actVelRE', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'nump', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', 'scalar'}));
addRequired(iP, 'pnames', ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addRequired(iP, 'plabels', ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addRequired(iP, 'pislog', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'pvalues', ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addRequired(iP, 'nperp', ...
    @(x) validateattributes(x, {'double'}, {'nonempty'}));
addRequired(iP, 'ncells', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', 'scalar'}));
addRequired(iP, 'actmode', ...
    @(x) validateattributes(x, {'double'}, {'nonempty', 'scalar'}));
addRequired(iP, 'loopmode', ...
    @(x) validateattributes(x, {'char' ,'string'}, {'nonempty'}));

  
% Add optional inputs to the Input Parser
%%% TODO

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigTypes', 'png', ... % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'CellsToPlot', 0:ncells/2-1, ...% cells to plot for latencies when not comparing across conditions
    @(x) validateattributes(x, {'numeric'}, {'vector', 'nonnegative', 'integer'}));
addParameter(iP, 'PCondNames', {}, ...  % names of parameters to be varied across conditions
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'PCondNames must be a cell array of strings or character arrays!'));

% Read from the Input Parser
parse(iP, infolders, outfolder, numActiveTC, numActiveRE, latencyTC, latencyRE, oscDurTC, oscDurRE, nSpikesTC, nSpikesRE, actDurTC, actDurRE, actVelTC, actVelRE, nump, pnames, plabels, pislog, pvalues, nperp, ncells, actmode, loopmode, varargin{:});
infolders = iP.Results.infolders;
outfolder = iP.Results.outfolder;
numActiveTC = iP.Results.numActiveTC;
numActiveRE = iP.Results.numActiveRE;
latencyTC = iP.Results.latencyTC;
latencyRE = iP.Results.latencyRE;
oscDurTC = iP.Results.oscDurTC;
oscDurRE = iP.Results.oscDurRE;
nSpikesTC = iP.Results.nSpikesTC;
nSpikesRE = iP.Results.nSpikesRE;
actDurTC = iP.Results.actDurTC;
actDurRE = iP.Results.actDurRE;
actVelTC = iP.Results.actVelTC;
actVelRE = iP.Results.actVelRE;
nump = iP.Results.nump;
pnames = iP.Results.pnames;
plabels = iP.Results.plabels;
pislog = iP.Results.pislog;
pvalues = iP.Results.pvalues;
nperp = iP.Results.nperp;
ncells = iP.Results.ncells;
actmode = iP.Results.actmode;
loopmode = iP.Results.loopmode;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
latency_cells_to_plot = iP.Results.CellsToPlot;
pcondnames = iP.Results.PCondNames;


%% Check if invariants are the same across conditions; if so, rename invariants
if across_conditions
    nump = check_and_collapse_identical_contents(nump, 'nump');
    pnames = check_and_collapse_identical_contents(pnames, 'pnames');
    plabels = check_and_collapse_identical_contents(plabels, 'plabels');
    pislog = check_and_collapse_identical_contents(pislog, 'pislog');
    ncells = check_and_collapse_identical_contents(ncells, 'ncells');
    actmode = check_and_collapse_identical_contents(actmode, 'actmode');
end

%% Create cell labels
cell_labels = cell(1, ncells);
for c = 1:ncells
    cell_labels{c} = ['Cell', num2str(c-1)];
end

%% Find varied parameters values for each condition
if across_conditions && isempty(pcondnames)
    pcondnames = cell(1, nump);
    if nump == 2            % currently only supports the case when 2 different parameters are looped and varied
        for p = 1:nump                  % for each parameter that is varied within a condition
            % Find parameter that is varied across conditions (i.e., the other parameter that is looped)
            q = 3 - p;                  % 2 -> 1; 1 -> 2 
            pcondnames{p} = pnames{q};  % the name of the parameter that is varied across conditions
        end
    else
        error('Please specify the names of the parameters varied across conditions!');
    end
end

%% Create condition labels
pcondvalues = cell(1, nump);
cond_labels = cell(1, nump);
cellcond_labels = cell(1, nump);
if across_conditions    
    for p = 1:nump          % for each parameter that is varied within a condition
        pcondvalues{p} = cell(1, ncond);
        cond_labels{p} = cell(1, ncond);
        cellcond_labels{p} = cell(1, ncond);
        for k = 1:ncond
            % Extract parameter values from sim_params file
            sim_params_files = dirr(fullfile(infolders{k}, ...
                ['sim_params_', pnames{p}, '*']), '.csv');    
                            % all sim_params files for the current looped parameter
            fid = fopen(fullfile(infolders{k}, sim_params_files(1).name));        % the first file suffices
            simfilecontent = textscan(fid, '%s %f %s', 'Delimiter', ',');
            paramnames = simfilecontent{1};
            params_val = simfilecontent{2};
            pcondvalues{p}{k} = params_val(find_in_strings(pcondnames{p}, paramnames));

            % Create condition label
            cond_labels{p}{k} = [pcondnames{p}, ' = ', num2str(pcondvalues{p}{k})];

            % Create cell-condition labels
            cellcond_labels{p}{k} = cell(1, ncells);
            for c = 1:ncells
                cellcond_labels{p}{k}{c} = ...
                    ['Cell', num2str(c-1), ', ', ...
                    pcondnames{p}, ' = ', num2str(pcondvalues{p}{k})];
            end            
        end
    end
end

%% Reorganize numActive, latency, oscDur 
%   so that each vector or matrix is for one parameter only
if across_conditions
    for k = 1:ncond
        numActiveTC{k} = reorganize_readout (numActiveTC{k}, nump, nperp{k});
	    numActiveRE{k} = reorganize_readout (numActiveRE{k}, nump, nperp{k});	
        latencyTC{k} = reorganize_readout (latencyTC{k}, nump, nperp{k});
		latencyRE{k} = reorganize_readout (latencyRE{k}, nump, nperp{k});
        oscDurTC{k} = reorganize_readout (oscDurTC{k}, nump, nperp{k});
        oscDurRE{k} = reorganize_readout (oscDurRE{k}, nump, nperp{k});
		nSpikesTC{k} = reorganize_readout (nSpikesTC{k}, nump, nperp{k});
		nSpikesRE{k} = reorganize_readout (nSpikesRE{k}, nump, nperp{k});
		actDurTC{k} = reorganize_readout (actDurTC{k}, nump, nperp{k});
		actDurRE{k} = reorganize_readout (actDurRE{k}, nump, nperp{k});
		actVelTC{k} = reorganize_readout (actVelTC{k}, nump, nperp{k});
		actVelRE{k} = reorganize_readout (actVelRE{k}, nump, nperp{k});
    end
else
    numActiveTC = reorganize_readout (numActiveTC, nump, nperp);
	numActiveRE = reorganize_readout (numActiveRE, nump, nperp);
    latencyTC = reorganize_readout (latencyTC, nump, nperp);
	latencyRE = reorganize_readout (latencyRE, nump, nperp);
    oscDurTC = reorganize_readout (oscDurTC, nump, nperp);
    oscDurRE = reorganize_readout (oscDurRE, nump, nperp);
    nSpikesTC = reorganize_readout (nSpikesTC, nump, nperp);
	nSpikesRE = reorganize_readout (nSpikesRE, nump, nperp);
	actDurTC = reorganize_readout (actDurTC, nump, nperp);
	actDurRE = reorganize_readout (actDurRE, nump, nperp);
	actVelTC = reorganize_readout (actVelTC, nump, nperp);
	actVelRE = reorganize_readout (actVelRE, nump, nperp);
end

%% Create tuning curves
for p = 1:nump
    % Define color map if many tuning curves are plotted
    if across_conditions
        cm = colormap(lines(ncond));    % color map used
    end
	
	if across_conditions
		if latency_cells_to_plot == -1
            % Plot the latencies of the center neuron by default
            cols_to_plot = ncells/2 + 1;
        else
            cols_to_plot = latency_cells_to_plot + 1;
        end
	else
		if latency_cells_to_plot == -1
            cols_to_plot = (0:(ncells/2-1)) + 1;
        else
            cols_to_plot = latency_cells_to_plot + 1;
        end
	end

    % Plot number of activated cells against parameter values
    ylimits = [0, ncells*2];
    if across_conditions
		figname = fullfile(outfolder, ['numActiveTC_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, numActiveTC{k}{p}, 1, ...
                        pislog(p), plabels{p}, numActiveTC_label, {cond_labels{p}{k}}, ...
                        [], ylimits, [], 'SingleColor', cm(k, :));        
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes);

		figname = fullfile(outfolder, ['numActiveRE_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, numActiveRE{k}{p}, 1, ...
                        pislog(p), plabels{p}, numActiveRE_label, {cond_labels{p}{k}}, ...
                        [], ylimits, [], 'SingleColor', cm(k, :));        
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes);
    else
		figname = fullfile(outfolder, ['numActiveTC_vs_', pnames{p}]);
        plot_tuning_curve(pvalues{p}, numActiveTC{p}, 1, ...
                    pislog(p), plabels{p}, numActiveTC_label, ...
                    {'suppress'}, [], ylimits, figname, ...
                    'FigTypes', figtypes);

		figname = fullfile(outfolder, ['numActiveRE_vs_', pnames{p}]);
        plot_tuning_curve(pvalues{p}, numActiveRE{p}, 1, ...
                    pislog(p), plabels{p}, numActiveRE_label, ...
                    {'suppress'}, [], ylimits, figname, ...
                    'FigTypes', figtypes);
        %% TODO: Pass 'figtypes' to plot_tuning_curve.m
    end

    % Plot latency to activation (seconds) of cells against parameter values
    if across_conditions
		figname = fullfile(outfolder, ['latencyTC_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1); % latencyTC
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, latencyTC{k}{p}, cols_to_plot, ...
                        pislog(p), plabels{p}, latencyTC_label, cellcond_labels{p}{k}, ...
                        [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)

		figname = fullfile(outfolder, ['latencyRE_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1); % latencyRE
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, latencyRE{k}{p}, cols_to_plot, ...
                        pislog(p), plabels{p}, latencyRE_label, cellcond_labels{p}{k}, ...
                        [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)
	else
		figname = fullfile(outfolder, ['latencyTC_vs_', pnames{p}]);
        yMaxTC = max(latencyTC{p}(~isinf(latencyTC{p}))) + 1;
        if isempty(yMaxTC) || yMaxTC < 1
            yMaxTC = 1;
        end
        ylimitsTC = [0, yMaxTC];
        plot_tuning_curve(pvalues{p}, latencyTC{p}, cols_to_plot, ... % latencyTC
                    pislog(p), plabels{p}, latencyTC_label, ...
                    cell_labels, [], ylimitsTC, figname, ...
                    'FigTypes', figtypes);
		
		figname = fullfile(outfolder, ['latencyRE_vs_', pnames{p}]);
        yMaxRE = max(latencyRE{p}(~isinf(latencyRE{p}))) + 1;
        if isempty(yMaxRE) || yMaxRE < 1
            yMaxRE = 1;
        end
        ylimitsRE = [0, yMaxRE];
        plot_tuning_curve(pvalues{p}, latencyRE{p}, cols_to_plot, ... % latencyRE
                    pislog(p), plabels{p}, latencyRE_label, ...
                    cell_labels, [], ylimitsRE, figname, ...
                    'FigTypes', figtypes);
	end
    % Plot oscillation duration against parameter values
    if across_conditions
		figname = fullfile(outfolder, ['oscDurTC_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, oscDurTC{k}{p}, 1, ...
                    pislog(p), plabels{p}, oscDurTC_label, {cond_labels{p}{k}}, ...
                    [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)

		figname = fullfile(outfolder, ['oscDurRE_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, oscDurRE{k}{p}, 1, ...
                    pislog(p), plabels{p}, oscDurRE_label, {cond_labels{p}{k}}, ...
                    [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)
    else
		figname = fullfile(outfolder, ['oscDurTC_vs_', pnames{p}]);
        yMax = max(oscDurTC{p})+0.1;
        if isempty(yMax)
            yMax = 1;
        end
        ylimits = [0, yMax];
        plot_tuning_curve(pvalues{p}, oscDurTC{p}, 1, ...
                pislog(p), plabels{p}, oscDurTC_label, ...
                {'suppress'}, [], ylimits, figname, ...
                'FigTypes', figtypes);

		figname = fullfile(outfolder, ['oscDurRE_vs_', pnames{p}]);
        yMax = max(oscDurRE{p})+0.1;
        if isempty(yMax)
            yMax = 1;
        end
        ylimits = [0, yMax];
        plot_tuning_curve(pvalues{p}, oscDurRE{p}, 1, ...
                pislog(p), plabels{p}, oscDurRE_label, ...
                {'suppress'}, [], ylimits, figname, ...
                'FigTypes', figtypes);
    end
	% Plot number of spikes against parameter values
	if across_conditions
		figname = fullfile(outfolder, ['nSpikesTC_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1); % nSpikesTC
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, nSpikesTC{k}{p}, cols_to_plot, ...
                        pislog(p), plabels{p}, nSpikesTC_label, cellcond_labels{p}{k}, ...
                        [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)

		figname = fullfile(outfolder, ['nSpikesRE_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1); % nSpikesRE
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, nSpikesRE{k}{p}, cols_to_plot, ...
                        pislog(p), plabels{p}, nSpikesRE_label, cellcond_labels{p}{k}, ...
                        [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)
	else
		figname = fullfile(outfolder, ['nSpikesTC_vs_', pnames{p}]);
        yMaxTC = max(nSpikesTC{p}(~isinf(nSpikesTC{p}))) + 1;
        if isempty(yMaxTC) || yMaxTC < 1
            yMaxTC = 1;
        end
        ylimitsTC = [0, yMaxTC];
        plot_tuning_curve(pvalues{p}, nSpikesTC{p}, cols_to_plot, ... % nSpikesTC
                    pislog(p), plabels{p}, nSpikesTC_label, ...
                    cell_labels, [], ylimitsTC, figname, ...
                    'FigTypes', figtypes);

		figname = fullfile(outfolder, ['nSpikesRE_vs_', pnames{p}]);
        yMaxRE = max(nSpikesRE{p}(~isinf(nSpikesRE{p}))) + 1;
        if isempty(yMaxRE) || yMaxRE < 1
            yMaxRE = 1;
        end
        ylimitsRE = [0, yMaxRE];
        plot_tuning_curve(pvalues{p}, nSpikesRE{p}, cols_to_plot, ... % nSpikesRE
                    pislog(p), plabels{p}, nSpikesRE_label, ...
                    cell_labels, [], ylimitsRE, figname, ...
                    'FigTypes', figtypes);
	end
	% Plot active duration (seconds) of cells against parameter values
	if across_conditions
		figname = fullfile(outfolder, ['actDurTC_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1); % actDurTC
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, actDurTC{k}{p}, cols_to_plot, ...
                        pislog(p), plabels{p}, actDurTC_label, cellcond_labels{p}{k}, ...
                        [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)

		figname = fullfile(outfolder, ['actDurRE_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1); % actDurRE
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, actDurRE{k}{p}, cols_to_plot, ...
                        pislog(p), plabels{p}, actDurRE_label, cellcond_labels{p}{k}, ...
                        [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)
    else
		figname = fullfile(outfolder, ['actDurTC_vs_', pnames{p}]);
        yMaxTC = max(actDurTC{p}(~isinf(actDurTC{p}))) + 1;
        if isempty(yMaxTC) || yMaxTC < 1
            yMaxTC = 1;
        end
        ylimitsTC = [0, yMaxTC];
        plot_tuning_curve(pvalues{p}, actDurTC{p}, cols_to_plot, ... % actDurTC
                    pislog(p), plabels{p}, actDurTC_label, ...
                    cell_labels, [], ylimitsTC, figname, ...
                    'FigTypes', figtypes);

		figname = fullfile(outfolder, ['actDurRE_vs_', pnames{p}]);
        yMaxRE = max(actDurRE{p}(~isinf(actDurRE{p}))) + 1;
        if isempty(yMaxRE) || yMaxRE < 1
            yMaxRE = 1;
        end
        ylimitsRE = [0, yMaxRE];
        plot_tuning_curve(pvalues{p}, actDurRE{p}, cols_to_plot, ... % actDurRE
                    pislog(p), plabels{p}, actDurRE_label, ...
                    cell_labels, [], ylimitsRE, figname, ...
                    'FigTypes', figtypes);
    end
    % Plot activation velocity (cells/second) against parameter values
    if across_conditions
		figname = fullfile(outfolder, ['actVelTC_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, actVelTC{k}{p}, 1, ...
                    pislog(p), plabels{p}, actVelTC_label, {cond_labels{p}{k}}, ...
                    [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)

		figname = fullfile(outfolder, ['actVelRE_vs_', pnames{p}]);
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, actVelRE{k}{p}, 1, ...
                    pislog(p), plabels{p}, actVelRE_label, {cond_labels{p}{k}}, ...
                    [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)
    else
		figname = fullfile(outfolder, ['actVelTC_vs_', pnames{p}]);
        yMax = max(actVelTC{p})+0.1;
        if isempty(yMax)
            yMax = 1;
        end
        ylimits = [0, yMax];
        plot_tuning_curve(pvalues{p}, actVelTC{p}, 1, ...
                pislog(p), plabels{p}, actVelTC_label, ...
                {'suppress'}, [], ylimits, figname, ...
                'FigTypes', figtypes);

		figname = fullfile(outfolder, ['actVelRE_vs_', pnames{p}]);
        yMax = max(actVelRE{p})+0.1;
        if isempty(yMax)
            yMax = 1;
        end
        ylimits = [0, yMax];
        plot_tuning_curve(pvalues{p}, actVelRE{p}, 1, ...
                pislog(p), plabels{p}, actVelRE_label, ...
                {'suppress'}, [], ylimits, figname, ...
                'FigTypes', figtypes);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readout_reorg = reorganize_readout (readout, nump, nperp)
%% Reorganize readout vector or matrix

%% Create class vector for the parameter classes
classvector = [];
for p = 1:nump
    classvector = [classvector; p * ones(nperp(p), 1)];
end

%% Reorganize readout vector into a cell array of vectors according to parameter class
readout_reorg = vec2cell(readout, classvector);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%