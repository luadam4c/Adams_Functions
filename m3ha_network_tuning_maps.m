function m3ha_network_tuning_maps (infolder, outfolder, numActiveTC, numActiveRE, ...
                latencyTC, latencyRE, oscDurTC, oscDurRE, ...
                nSpikesTC, nSpikesRE, actDurTC, actDurRE, ...
                actVelTC, actVelRE, nump, pnames, plabels, ...
                pislog, pvalues, nperp, ncells, actmode, loopmode, varargin)
%% Shows a tuning map for numActive, latency, oscDur, for each parameter changed
% USAGE: m3ha_network_tuning_maps (infolder, outfolder, numActiveTC, numActiveRE, ...
%               latencyTC, latencyRE, oscDurTC, oscDurRE, ...
%               nSpikesTC, nSpikesRE, actDurTC, actDurRE, ...
%               actVelTC, actVelRE, nump, pnames, plabels, ...
%               pislog, pvalues, nperp, ncells, actmode, loopmode, varargin)
% Arguments:
%       infolder    - the directory that contains the data to plot
%       outfolder   - (opt) the directory to place the output plots
%           TODO: other optional arguments
%       varargin    - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - 'CellToPlot': cell to plot for the latency tuning map
%                   default == ncells/2
%
% Requires:
%        cd/extract_looped_params.m
%        cd/isfigtype.m
%        cd/m3ha_network_raster_plot.m
%        cd/plot_tuning_map.m
%        cd/save_all_figtypes.m
%        cd/vec2array.m
%        infolder/*loopedparams.mat
%
% Used by:
%        cd/m3ha_network_raster_plot.m
%
% 2017-10-23 Modified from /RTCl/m3ha_network_tuning_maps.m
% 2018-01-26 BT - Differentiate between TC and RE
% 2018-01-31 BT - Updated with numActive, oscDur, and actVel TC/RE

%% Set readout labels
numActiveTC_label = 'Number of activated TC cells';
latencyTC_label = 'TC Latency to activation (seconds)';
oscDurTC_label = 'TC Oscillation duration (seconds)';
nSpikesTC_label = 'TC Number of spikes';
actDurTC_label = 'TC Active duration (seconds)';
actVelTC_label = 'TC Active velocity (cells/second)';
numActiveRE_label = 'Number of activated RE cells';
latencyRE_label = 'RE Latency to activation (seconds)';
oscDurRE_label = 'RE Oscillation duration (seconds)';
nSpikesRE_label = 'RE Number of spikes';
actDurRE_label = 'RE Active duration (seconds)';
actVelRE_label = 'RE Active velocity (cells/second)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
addpath(fullfile(functionsdirectory, '/Adams_Functions/'));        % for isfigtype.m, extract_looped_params.m,
                                    %     vec2array.m
                                    %    & plot_tuning_map.m
%% Check arguments
%%% TODO

%% Set outfolder and create if doesn't exist
if nargin < 2 || isempty(outfolder)
    outfolder = infolder;
end
if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
end

%% Get numActive, latency, oscDur if not provided
if nargin < 3 || isempty(numActiveTC) || isempty(numActiveRE) || isempty(latencyTC) || isempty(latencyRE) || isempty(oscDurTC) || isempty(oscDurRE) % update with input parser
	[~, ~, numActive, latencyTC, latencyRE, oscDur, nSpikesTC, nSpikesRE, actDurTC, actDurRE, actVel] = ...
            m3ha_network_raster_plot(infolder, 'OutFolder', outfolder, ...
                'RenewParpool', 0, 'PlotSpikes', 0, 'PlotTuning', 0);
end

%% Get looped parameters info if not provided
if nargin < 6 || isempty(nump) || isempty(pnames) || isempty(plabels) ...
        || isempty(pislog) || isempty(pvalues) || isempty(nperp) ...
        || isempty(ncells) || isempty(actmode) || isempty(loopmode)
    [nump, pnames, plabels, pislog, pvalues, nperp, ~, ~, ncells, actmode, loopmode] = extract_looped_params(infolder);
end

%% Use Input Parser for parameter-value pairs
iP = inputParser;
addParameter(iP, 'CellToPlot', -1);         % cell to plot for the latency tuning map
addParameter(iP, 'FigTypes', 'png', ...     % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
parse(iP, varargin{:});
latency_cell_to_plot = iP.Results.CellToPlot;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Find latency of just the cell to plot
if latency_cell_to_plot == -1
    % Plot the latencies of the center TC neuron by default
    latency_cell_to_plot = ncells/2;
end
latencyTC_this = latencyTC(:, latency_cell_to_plot + 1);
latencyRE_this = latencyRE(:, latency_cell_to_plot + 1);

nSpikesTC_this = nSpikesTC(:, latency_cell_to_plot + 1);
nSpikesRE_this = nSpikesRE(:, latency_cell_to_plot + 1);
actDurTC_this = actDurTC(:, latency_cell_to_plot + 1);
actDurRE_this = actDurRE(:, latency_cell_to_plot + 1);

%% Reorganize numActive, latency, oscDur so that it is an array with nump dimensions
numActiveTC = vec2array(numActiveTC, nperp);
numActiveRE = vec2array(numActiveRE, nperp);
latencyTC_this = vec2array(latencyTC_this, nperp);
latencyRE_this = vec2array(latencyRE_this, nperp);
oscDurTC = vec2array(oscDurTC, nperp);
oscDurRE = vec2array(oscDurRE, nperp);
nSpikesTC_this = vec2array(nSpikesTC_this, nperp);
nSpikesRE_this = vec2array(nSpikesRE_this, nperp);
actDurTC_this = vec2array(actDurTC_this, nperp);
actDurRE_this = vec2array(actDurRE_this, nperp);
actVelTC = vec2array(actVelTC, nperp);
actVelRE = vec2array(actVelRE, nperp);

%% Create heat maps if there are two parameters
if nump == 2
    % Plot number of activated cells against parameter values
    figname = fullfile(outfolder, ['numActiveTC_vs_', pnames{1}, '_', pnames{2}]);
    climits = [0, ncells];
    plot_tuning_map (pvalues, numActiveTC, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', numActiveTC_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    figname = fullfile(outfolder, ['numActiveRE_vs_', pnames{1}, '_', pnames{2}]);
    climits = [0, ncells];
    plot_tuning_map (pvalues, numActiveRE, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', numActiveRE_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)
    % Plot latency to activation (seconds) of cells against parameter values
    figname = fullfile(outfolder, ['latencyTC_vs_', pnames{1}, '_', pnames{2}]);
    maxLatency = max(1, max(max(latencyTC(~isinf(latencyTC)))));
    if ~isnan(maxLatency)
        yMax = maxLatency + 1;
    else
        yMax = 1;
    end	
    climits = [0, yMax];
    latencyTC_label_this = [latencyTC_label, ' for Cell#', num2str(latency_cell_to_plot)];
    plot_tuning_map (pvalues, latencyTC_this, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', latencyTC_label_this, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    figname = fullfile(outfolder, ['latencyRE_vs_', pnames{1}, '_', pnames{2}]);
    maxLatency = max(1, max(max(latencyRE(~isinf(latencyRE)))));
    if ~isnan(maxLatency)
        yMax = maxLatency + 1;
    else
        yMax = 1;
    end
    climits = [0, yMax];
    latencyRE_label_this = [latencyRE_label, ' for Cell#', num2str(latency_cell_to_plot)];
    plot_tuning_map (pvalues, latencyRE_this, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', latencyRE_label_this, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    % If not stimulating everywhere, plot oscillation duration (seconds) against parameter values
    figname = fullfile(outfolder, ['oscDurTC_vs_', pnames{1}, '_', pnames{2}]);
    maxOscDur = max(max(abs(oscDurTC)));
    if ~isnan(maxOscDur)
        yMax = maxOscDur + 0.1;
    else
        yMax = 1;
    end
    climits = [0, yMax];
    plot_tuning_map (pvalues, oscDurTC, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', oscDurTC_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    figname = fullfile(outfolder, ['oscDurRE_vs_', pnames{1}, '_', pnames{2}]);
    maxOscDur = max(max(abs(oscDurRE)));
    if ~isnan(maxOscDur)
        yMax = maxOscDur + 0.1;
    else
        yMax = 1;
    end
    climits = [0, yMax];
    plot_tuning_map (pvalues, oscDurRE, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', oscDurRE_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

	% Plot number of spikes of cells against parameter values
    figname = fullfile(outfolder, ['nSpikesTC_vs_', pnames{1}, '_', pnames{2}]);
    maxNSpikes = max(1, max(max(nSpikesTC(~isinf(nSpikesTC)))));
    if ~isnan(maxNSpikes)
        yMax = maxNSpikes + 1;
    else
        yMax = 1;
    end	
    climits = [0, yMax];
    plot_tuning_map (pvalues, nSpikesTC_this, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', nSpikesTC_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    figname = fullfile(outfolder, ['nSpikesRE_vs_', pnames{1}, '_', pnames{2}]);
    maxNSpikes = max(1, max(max(nSpikesRE(~isinf(nSpikesRE)))));
    if ~isnan(maxNSpikes)
        yMax = maxNSpikes + 1;
    else
        yMax = 1;
    end
    climits = [0, yMax];
    plot_tuning_map (pvalues, nSpikesRE_this, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', nSpikesRE_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

	% Plot active duration (seconds) of cells against parameter values
    figname = fullfile(outfolder, ['actDurTC_vs_', pnames{1}, '_', pnames{2}]);
    maxNSpikes = max(1, max(max(actDurTC(~isinf(actDurTC)))));
    if ~isnan(maxNSpikes)
        yMax = maxNSpikes + 1;
    else
        yMax = 1;
    end	
    climits = [0, yMax];
    plot_tuning_map (pvalues, actDurTC_this, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', actDurTC_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    figname = fullfile(outfolder, ['actDurRE_vs_', pnames{1}, '_', pnames{2}]);
    maxNSpikes = max(1, max(max(actDurRE(~isinf(actDurRE)))));
    if ~isnan(maxNSpikes)
        yMax = maxNSpikes + 1;
    else
        yMax = 1;
    end
    climits = [0, yMax];
    plot_tuning_map (pvalues, actDurRE_this, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', actDurRE_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    % Plot activation velocity (seconds/cell) against parameter values
    figname = fullfile(outfolder, ['actVelTC_vs_', pnames{1}, '_', pnames{2}]);
    maxActVel = max(max(abs(actVelTC)));
    if ~isnan(maxActVel)
        yMax = maxActVel + 0.1;
    else
        yMax = 1;
    end
    climits = [0, yMax];
    plot_tuning_map (pvalues, actVelTC, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', actVelTC_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    figname = fullfile(outfolder, ['actVelRE_vs_', pnames{1}, '_', pnames{2}]);
    maxActVel = max(max(abs(actVelRE)));
    if ~isnan(maxActVel)
        yMax = maxActVel + 0.1;
    else
        yMax = 1;
    end
    climits = [0, yMax];
    plot_tuning_map (pvalues, actVelRE, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', actVelRE_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
