function [data] = m3ha_network_single_neuron (inFolder, varargin)
%% Shows single neuron traces for different neurons or for different properties in the same neuron
% USAGE: [data] = m3ha_network_single_neuron (inFolder, varargin)
% Arguments:
%   inFolder    - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%            must be a character array
%   varargin    - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%            could be anything recognised by the built-in saveas() function
%            (see isfigtype.m under Adams_Functions)
%            default == 'png'
%            - 'OutFolder': the name of the directory that the plots will be placed
%            must be a directory
%            default: same as inFolder
%            - 'MaxNumWorkers': maximum number of workers for running NEURON 
%            must a positive integer
%            default: 20
%            - 'RenewParpool': whether to renew parpool every batch to release memory
%            must be logical 1 (true) or 0 (false)
%            default: true, but alway false if plotspikes == false
%             - 'CellsToPlot': the ID #s for cells whose voltage & chloride concentration traces are to be plotted
%            must be a numeric array with elements that are integers between 0 and nCells
%            default: [act, actLeft1, actLeft2, far], whose values are saved in sim_params.csv
%            - 'PropertiesToPlot': property #s of special neuron to record to be plotted 
%            maximum range: 1~12, must be consistent with m3ha_net.hoc
%            must be a numeric array with elements that are integers between 0 and 12
%            default: 1:12
%            legend for RTCl:    TODO: for m3ha
%                       1 - voltage (mV) trace
%                       2 - sodium current (mA/cm2) trace
%                       3 - potassium current (mA/cm2) trace
%                       4 - calcium current (mA/cm2) trace 
%                       5 - calcium concentration (mM) trace
%                       6 - GABA-A chloride current (nA) trace
%                       7 - GABA-A bicarbonate current (nA) trace
%                       8 - chloride current (mA/cm2) trace
%                       9 - chloride concentration (mM) trace
%                       10 - chloride concentration (mM) in inner annuli trace
%                       11 - chloride reversal potential trace
%                       12 - GABA-A reversal potential trace
%            NOTE: must be consistent with propLabels & m3ha_net.hoc
%
% Requires:
%       cd/all_files.m
%       cd/extract_fileparts.m
%       cd/find_in_strings.m
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%       inFolder/*.singv OR inFolder/*.singcli OR inFolder/*.singsp
%       inFolder/['sim_params_', pstring, '.csv'] for all the possible parameter strings
%       /home/Matlab/Downloaded_Functions/subaxis.m
%
% Used by:
%       cd/m3ha_launch.m

% File History:
% 2017-10-23 Modified from /RTCl/single_neuron.m
% 2017-10-31 Replaced REuseca with useHH
% 2017-11-03 Changed ylim for voltage [-100, 60] -> [-120, 60]
% 2018-04-17 Fixed legend() to conform with R2017a
% 2018-04-27 Now plots spikes on voltage traces of special neuron plots
% 2019-11-04 Updated code to use all_files.m, extract_fileparts.m
% TODO: Take useHH as an optional argument and change propLabelsTC accordingly
%

%% Hard-coded parameters
% Set parameters
nZooms = 4;                 % number of different time intervals to plot for each data
nPlotsPerFile = nZooms + 1; % number of plots per file

% Set property labels
%   Note: must be consistent with m3ha_net.hoc
propLabelsRT = {'v (mV)', 'ina (mA/cm2)', 'ik (mA/cm2)', ...
                'ica (mA/cm2)', 'iAMPA (nA)', 'iGABA (nA)', ...
                'cai (mM)', 'cli (mM)', 'Gicl (nA)', 'Gihco3 (nA)', ...
                'icl (mA/cm2)', 'cli1 (mM)', ...
                'ecl (mV)', 'eGABA (mV)'};
propLabelsTC = {'v (mV)', 'ina (mA/cm2)', 'ik (mA/cm2)', ...
                'ica (mA/cm2)', 'iGABA_A (nA)', 'iGABA_B (nA)', ...
                'cai (mM)', 'gGABA_B (uS)', ...
                'm_{T,dend2}', 'm_{\infty,T,dend2}', ...
                'h_{T,dend2}', 'h_{\infty,T,dend2}'};
% propLabelsTC = {'v (mV)', 'ina (mA/cm2)', 'ik (mA/cm2)', ...
%                 'ica (mA/cm2)', 'iGABA_A (nA)', 'iGABA_B (nA)', ...
%                 'cai (mM)', 'gGABA_B (uS)'};
% propLabelsTC = {'v (mV)', 'inRefractory', 'ik (mA/cm2)', ...
%                 'ica (mA/cm2)', 'iGABA_A (nA)', 'iGABA_B (nA)', ...
%                 'cai (mM)', 'gGABA_B (uS)'};
% propLabelsTC = {'v (mV)', 'inSlopeWatching', 'ik (mA/cm2)', ...
%                 'ica (mA/cm2)', 'iGABA_A (nA)', 'iGABA_B (nA)', ...
%                 'cai (mM)', 'gGABA_B (uS)'};

% Set figure name suffixes
figSuffixVoltage = '_selected_soma_voltage.png';
figSuffixCli = '_selected_soma_cli.png';
figSuffixSp = '_alltraces.png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'inFolder', @isdir);    % the name of the directory containing the .syn files

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigTypes', 'png', ... % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'OutFolder', '@inFolder', @isdir);     % the name of the directory that the plots will be placed
addParameter(iP, 'MaxNumWorkers', 20, ...               % maximum number of workers for running NEURON
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'RenewParpool', true, ...              % whether to renew parpool every batch to release memory
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CellsToPlot', [], ...                 % the ID #s for cells to be plotted
    @(x) validateattributes(x, {'numeric'}, {'vector', 'nonnegative', 'integer', '>=', 0}));
addParameter(iP, 'PropertiesToPlot', 1:1:12, ...        % property #s of special neurons to be plotted
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer', ...
                            '>', 0, '<', 13}));
addParameter(iP, 'PropertiesToPlotRT', [], ...          % property #s of special RT neurons to be plotted
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer', ...
                            '>', 0, '<', 15}));
addParameter(iP, 'PropertiesToPlotTC', [], ...          % property #s of special TC neurons to be plotted
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer', ...
                            '>', 0, '<', 13}));

% Read from the Input Parser
parse(iP, inFolder, varargin{:});
[~, figtypes]    = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
outFolder        = iP.Results.OutFolder;
maxNumWorkers    = iP.Results.MaxNumWorkers;
renewParpool     = iP.Results.RenewParpool;
cellsToPlot      = iP.Results.CellsToPlot;
propertiesToPlot = iP.Results.PropertiesToPlot;
propertiesToPlotRT = iP.Results.PropertiesToPlotRT;
propertiesToPlotTC = iP.Results.PropertiesToPlotTC;

% Change default arguments if necessary
if strcmp(outFolder, '@inFolder')
    outFolder = inFolder;
end

%% Preparation
% Set default properties to plot
if isempty(propertiesToPlotRT)
    propertiesToPlotRT = propertiesToPlot;
end
if isempty(propertiesToPlotTC)
    propertiesToPlotTC = propertiesToPlot;
end

%% Find all .singv, .singcli or .singsp files
% Find files
[~, singPaths] = all_files('Directory', inFolder, 'Regexp', '.*.sing.*');

% Count the number of files
nFiles = numel(singPaths);

%% Set up plots
nPlots = nFiles * nPlotsPerFile;
fileNames = cell(nPlots, 1);    % stores file name for each plot
fileTypes = cell(nPlots, 1);    % stores file type for each plot
fullFigNames = cell(nPlots, 1); % stores full figure names for each plot
nCells = zeros(nPlots, 1);      % stores nCells for each plot
tStarts = zeros(nPlots, 1);     % stores tStart for each plot
tStops = zeros(nPlots, 1);      % stores tStop for each plot
stimStarts = zeros(nPlots, 1); % stores stimStart for each plot
stimDurs = zeros(nPlots, 1);   % stores stimDur for each plot
stimFreqs = zeros(nPlots, 1);  % stores stimFreq for each plot
useHH = zeros(nPlots, 1);       % stores useHH for each plot
act = zeros(nPlots, 1);         % stores act for each plot
actLeft1 = zeros(nPlots, 1);   % stores actLeft1 for each plot
actLeft2 = zeros(nPlots, 1);   % stores actLeft2 for each plot
far = zeros(nPlots, 1);         % stores far for each plot
cellIDsToPlot = cell(nPlots, 1);     % stores default neuron ID #s to plot for each plot
propLabels = cell(nFiles, 1);   % stores property labels for each plot
propertiesToPlotAll = cell(nFiles, 1);
nPlotsSet = 0;
for iFile = 1:nFiles
    % Get the current file path
    singPathThis = singPaths{iFile};

    % Extract the file base and extension
    fileName = extract_fileparts(singPathThis, 'name');
    fileBase = extract_fileparts(singPathThis, 'base');
    fileExt = extract_fileparts(singPathThis, 'extension');

    % Set filetype according to filename
    switch fileExt
        case '.singv'
            fileType = 'v';
        case '.singcli'
            fileType = 'cli';
        case '.singsp'
            fileType = 'sp';
        otherwise
            error('fileExt unrecognised!');
    end

    % Construct path to the simulation parameters spreadsheet
    simParamsFileName = ['sim_params_', extractAfter(fileBase, '_'), '.csv'];
    simParamsPath = fullfile(inFolder, simParamsFileName);

    % Read the simulation parameters table
    simParamsTable = readtable(simParamsPath, 'ReadRowNames', true);

    % Set property labels according to file name
    switch fileBase(1:2)
    case 'RE'
        propLabels{iFile} = propLabelsRT;
        propertiesToPlotAll{iFile} = propertiesToPlotRT;
    case 'TC'
        propLabels{iFile} = propLabelsTC;
        propertiesToPlotAll{iFile} = propertiesToPlotTC;
    otherwise
        error('File name must include ''RE'' or ''TC''!');
    end

    % Find indices of plots corresponding to this file
    indPlots = nPlotsSet + (1:nPlotsPerFile);

    % Extract parameters from the simulation parameters table
    [nCellsThis, tStartThis, tStopThis, ...
        stimStartThis, stimDurThis, stimFreqThis, ...
        useHHThis, actThis, actLeft1This, ...
        actLeft2This, farThis] = ...
            argfun(@(x) simParamsTable{x, 'Value'}, ...
                    'nCells', 'tStart', 'tStop', ...
                    'stimStart', 'stimDur', 'stimFreq', ...
                    'useHH', 'act', 'actLeft1', ...
                    'actLeft2', 'far');

    % Store file names and types
    fileNames(indPlots) = repmat({fileName}, nPlotsPerFile, 1);
    fileTypes(indPlots) = repmat({fileType}, nPlotsPerFile, 1);

    % Store parameters
    [nCells(indPlots), stimStarts(indPlots), ...
        stimDurs(indPlots), stimFreqs(indPlots), ...
        useHH(indPlots), act(indPlots), actLeft1(indPlots), ...
        actLeft2(indPlots), far(indPlots)] = ...
            argfun(@(x) repmat(x, nPlotsPerFile, 1), ...
                    nCellsThis, stimStartThis, ...
                    stimDurThis, stimFreqThis, ...
                    useHHThis, actThis, actLeft1This, ...
                    actLeft2This, farThis);

    % Set default ID #s for neurons whose voltage is to be plotted
    if ~isempty(cellsToPlot)
        cellIDsToPlot(indPlots) = repmat({cellsToPlot}, nPlotsPerFile, 1);
    else
        cellIDsToPlot(indPlots) = ...
            arrayfun(@(x) [act(x), actLeft1(x), actLeft2(x), far(x)], ...
                    indPlots, 'UniformOutput', false);
    end

    % Set general figure names according to fileType
    switch fileType
        case 'v'
            figSuffix = figSuffixVoltage;
        case 'cli'
            figSuffix = figSuffixCli;
        case 'sp'
            figSuffix = figSuffixSp;
    end
    figName = [fileBase, figSuffix];
    figPath = fullfile(outFolder, figName);

    % Create full figure names with modifications
    fullFigNames{nPlotsSet+1} = figPath;
    fullFigNames{nPlotsSet+2} = replace(figPath, '.png', '_zoom1.png');
    fullFigNames{nPlotsSet+3} = replace(figPath, '.png', '_zoom2.png');
    fullFigNames{nPlotsSet+4} = replace(figPath, '.png', '_zoom3.png');
    fullFigNames{nPlotsSet+5} = replace(figPath, 'selected', 'heatmap');

    % Create time limits for different time intervals to plot
    %   TODO: Change for m3ha
    tStarts(nPlotsSet+1) = tStartThis;                                          % 0 ms
    tStops(nPlotsSet+1) = tStopThis;                                            % 30000 ms
    tStarts(nPlotsSet+2) = max(stimStartThis*2/3, tStartThis);                  % 200 ms
    tStops(nPlotsSet+2) = min(max(4000, stimStartThis*40/3), tStopThis);        % 30000 ms
    tStarts(nPlotsSet+3) = max(stimStartThis - 100, tStartThis);                % 2900 ms
    tStops(nPlotsSet+3) = min(stimStartThis + 400, tStopThis);                  % 4000 ms
    tStarts(nPlotsSet+4) = max(stimStartThis + stimDurThis, tStartThis);        % 500 ms
    tStops(nPlotsSet+4) = min(stimStartThis + stimDurThis + 1000, tStopThis);   % 1000 ms
%{
    tStarts(nPlotsSet+4) = max(stimStartThis + stimDurThis - 100, tStartThis);  % 500 ms
    tStops(nPlotsSet+4) = min(stimStartThis + stimDurThis + 400, tStopThis);    % 1000 ms
%}
    tStarts(nPlotsSet+5) = tStartThis;                                              % 0 ms
    tStops(nPlotsSet+5) = tStopThis;                                                % 233000 ms

    % Update nPlotsSet
    nPlotsSet = nPlotsPerFile .* iFile;
end

%% Create plots
data = cell(nPlots, 1);     % some elements will be empty but the indexing is necessary for parfor
ct = 0;                     % counts number of trials completed
poolObj = gcp('nocreate');  % get current parallel pool object without creating a new one
if isempty(poolObj)
    poolObj = parpool;      % create a default parallel pool object
    oldNumWorkers = poolObj.NumWorkers;         % number of workers in the default parallel pool object
else
    oldNumWorkers = poolObj.NumWorkers;         % number of workers in the current parallel pool object
end
numWorkers = min(oldNumWorkers, maxNumWorkers); % number of workers to use for running NEURON
if renewParpool
    delete(poolObj);        % delete the parallel pool object to release memory
end
while ct < nPlots           % while not trials are completed yet
    first = ct + 1;         % first trial in this batch
    if renewParpool && ct + numWorkers <= nPlots% if memory is to be released
        last = ct + numWorkers;                 % limit the batch to numWorkers
    else
        last = nPlots;
    end
    if renewParpool
        poolObj = parpool('local', numWorkers); % recreate a parallel pool object 
                            % using fewer workers to prevent running out of memory
    end

    parfor k = first:last
        iFile = ceil(k/nPlotsPerFile);
        if strcmp(fileTypes{k}, 'v') || strcmp(fileTypes{k}, 'cli')
            % Plot voltage or chloride concentration traces
            if mod(k, nPlotsPerFile) == 0
                plot_heat_map(fileTypes{k}, tStarts(k), tStops(k), fileNames{k}, fullFigNames{k}, ...
                        stimStarts(k), stimDurs(k), stimFreqs(k), inFolder, figtypes);
            elseif mod(k, nPlotsPerFile) == 1
                % data{k} = ...
                    plot_single_neuron_data(fileTypes{k}, nCells(k), useHH(k), ...
                        tStarts(k), tStops(k), fileNames{k}, fullFigNames{k}, ...
                        cellIDsToPlot{k}, stimStarts(k), stimDurs(k), stimFreqs(k), ...
                        inFolder, figtypes, propLabels{iFile});
            else
                plot_single_neuron_data(fileTypes{k}, nCells(k), useHH(k), ...
                    tStarts(k), tStops(k), fileNames{k}, fullFigNames{k}, ...
                    cellIDsToPlot{k}, stimStarts(k), stimDurs(k), stimFreqs(k), ...
                    inFolder, figtypes, propLabels{iFile});
            end
        elseif strcmp(fileTypes{k}, 'sp')
            % Plot other properties traces for special neurons
            if mod(k, nPlotsPerFile) == 0
                % No heat map; do nothing
            elseif mod(k, nPlotsPerFile) == 1
                % data{k} = ...
                    plot_single_neuron_data(fileTypes{k}, nCells(k), useHH(k), ...
                        tStarts(k), tStops(k), fileNames{k}, fullFigNames{k}, ...
                        propertiesToPlotAll{iFile}, stimStarts(k), stimDurs(k), stimFreqs(k), ...
                        inFolder, figtypes, propLabels{iFile});
            else
                plot_single_neuron_data(fileTypes{k}, nCells(k), useHH(k), ...
                    tStarts(k), tStops(k), fileNames{k}, fullFigNames{k}, ...
                    propertiesToPlotAll{iFile}, stimStarts(k), stimDurs(k), stimFreqs(k), ...
                    inFolder, figtypes, propLabels{iFile});
            end
        end

        close all;
    end
    if renewParpool
        delete(poolObj);    % delete the parallel pool object to release memory
    end
    ct = last;              % update number of trials completed
end
if renewParpool
    poolObj = parpool('local', oldNumWorkers);    % recreate a parallel pool object using the previous number of workers
end

%% Remove empty elements from data
data = data(~cellfun(@isempty, data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = plot_heat_map(filetype, tStart, tStop, filename, figName, stimStart, stimDur, stimFreq, inFolder, figtypes)
%% Plot heat map

% Load single neuron data
data = load(fullfile(inFolder, filename));

% Check (ID #s of neurons) to plot
nCells = size(data, 2) - 1;        % total number of columns in the data minus the time vector
if nCells < 1
    fprintf('Warning: No neurons to plot for this file!\n');
    return;
end

% Find range of data values to plot and load spike data
[climits, spikes] = get_aux(filetype, filename, inFolder);

% Get the spike cell numbers if available
if ~isempty(spikes)
    spikecelln = spikes(:, 1);
end

% Change units of time axis from ms to s if the total time is > 10 seconds
if ~isempty(spikes)
    [timevec, timelabel, xlim1, xlim2, stimStartPlot, stimDurPlot, spiketimes] ...
        = set_time_units(data, tStart, tStop, stimStart, stimDur, spikes);
else
    [timevec, timelabel, xlim1, xlim2, stimStartPlot, stimDurPlot, ~] ...
        = set_time_units(data, tStart, tStop, stimStart, stimDur);
end

% Find maximum and minimum time to plot
xmin = min([xlim1, min(timevec)]);
xmax = max([xlim2, max(timevec)]);

% Create plot
%h = figure(floor(rand()*10^4+, 1));
h = figure(10000);
clf(h);
hold on;

% Create heat map
imagesc([min(timevec), max(timevec)], [nCells-1, 0], flipud(data(:, 2:end)'));
set(gca, 'CLim', climits);
%HeatMap(flipud(data(:, 2:end)'), 'ColumnLabels', timevec, 'RowLabels', 0:nCells-1);        % Doesn't seem to work
%heatmap(timevec, 0:nCells-1, flipud(data(:, 2:end)'));     % Not available until R2017a
if ~isempty(spikes)
    plot(spiketimes, spikecelln, 'r.', 'MarkerSize', 1);    % plot spikes
end
line([stimStartPlot, stimStartPlot], [-1, nCells], ...
    'Color', 'r', 'LineStyle', '--');   % line for stimulation on
text(stimStartPlot + 0.5, nCells*0.95, ['Stim ON: ', num2str(stimFreq), ' Hz'], ...
    'Color', 'r');                      % text for stimulation on
line([stimStartPlot + stimDurPlot, stimStartPlot + stimDurPlot], ...
    [-1,  nCells], 'Color', 'r', 'LineStyle', '--');    % line for stimulation off
text(stimStartPlot + stimDurPlot + 0.5, nCells*0.95, 'Stim OFF', ...
    'Color', 'r');                      % text for stimulation off
xlim([xmin, xmax]);                     % time range to plot
ylim([-1, nCells]);                     % cell ID runs from 0 to nCells-1
xlabel(timelabel);
ylabel('Neuron number');
colorbar;
if strcmp(filetype, 'v')
    title(['Somatic voltage (mV) for ', replace(filename, '_', '\_')]);
elseif strcmp(filetype, 'cli')
    title(['Chloride concentration (mM) for ', replace(filename, '_', '\_')]);
end

% Save figure
save_all_figtypes(h, figName, figtypes);
%close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = plot_single_neuron_data(filetype, nCells, useHH, tStart, ...
                                        tStop, filename, figName, ToPl, ...
                                        stimStart, stimDur, stimFreq, ...
                                        inFolder, figtypes, propLabels)
%% Plot single neuron data

%% Extract info from arguments
nsubplots = length(ToPl);                % number of column vectors to plot

%% Create legend labels
if strcmp(filetype, 'v') || strcmp(filetype, 'cli')
    % Check ID numbers
    if max(ToPl) > nCells
        error('IDs are out of range!');
    end

    % Create ID labels
    labels = cell(1, nCells);    % ID labels
    for id = 0:nCells-1
        labels{id+1} = sprintf('Cell #%d', id);
    end
elseif strcmp(filetype, 'sp')
    % Set labels for properties of the special neuron to be plotted, 
    %       must be consistent with m3ha_net.hoc
    labels = propLabels;
end

% Load data
data = load(fullfile(inFolder, filename));

% Check (ID #s of neurons) or (properties of special neuron) to plot
ncols = size(data, 2) - 1;        % total number of columns in the data minus the time vector
if ncols < 1
    fprintf('Warning: No neurons or properties to plot for this file!\n');
    return;
end
for k = 1:nsubplots
    if strcmp(filetype, 'v') || strcmp(filetype, 'cli')
        if ToPl(k) < 0 || ToPl(k) > ncols - 1
            fprintf('Warning: ToPl(%d) is out of range; plotting first neuron or property instead\n', k);
            ToPl(k) = 0;
        end
    elseif strcmp(filetype, 'sp')
        if ToPl(k) < 1 || ToPl(k) > ncols
            fprintf('Warning: ToPl(%d) is out of range; plotting first neuron or property instead\n', k);
            ToPl(k) = 0;
        end
    end
end

% Get spikes for special neuron plots
% TODO: Do this for voltage plots too
if strcmp(filetype, 'sp')
    % Load spike data for this condition
    [~, spikes] = get_aux(filetype, filename, inFolder);

    % If spikes are available, get the spikes for this special neuron
    if ~isempty(spikes)
        % Get the base of the file name
        [~, fileBase, ~] = fileparts(filename);

        % Find the cell ID of interest
        temp1 = strsplit(fileBase, ']');
        temp2 = strsplit(temp1{1}, '[');
        cellID = str2double(temp2{2});

        % Get all spike cell IDs
        spikecelln = spikes(:, 1);

        % Get all indices corresponding to the cell of interest
        indThis = find(spikecelln == cellID);
    else
        indThis = [];
    end
end

% Deal with times
if strcmp(filetype, 'sp') && ~isempty(indThis)
    % Change units of time axis from ms to s if the total time is > 10 seconds
    [timevec, timelabel, xlim1, xlim2, stimStartPlot, stimDurPlot, spiketimes] ...
        = set_time_units(data, tStart, tStop, stimStart, stimDur, spikes);

    % Get the spike times for this cell
    spiketimesThis = spiketimes(indThis);

    % Find the corresponding voltage values for each neuron
    voltageThis = zeros(size(spiketimesThis));
    for iSpike = 1:length(spiketimesThis)
        % Find the last index of the time vector before spike time
        idxTime = find(timevec <= spiketimesThis(iSpike), 1, 'last');
        if isempty(idxTime)
            idxTime = 1;
        end

        % Find the corresponding voltage value
        voltageThis(iSpike) = data(idxTime, ToPl(1) + 1);
    end
else
    % Change units of time axis from ms to s if the total time is > 10 seconds
    [timevec, timelabel, xlim1, xlim2, stimStartPlot, stimDurPlot, ~] ...
        = set_time_units(data, tStart, tStop, stimStart, stimDur);
end

% Create figure
%h = figure(floor(rand()*10^4+, 1));
h = figure(10000);
clf(h);

% Plot voltage trace for each neuron with iD # in ToPl
for iSubPlot = 1:nsubplots
    % Generate subplot
    subaxis(nsubplots, 1, iSubPlot, 'SpacingVert', 0.015)
    hold on;
    if strcmp(filetype, 'v') || strcmp(filetype, 'cli')
        % Trace for neuron #i is in the i+2nd column
        %   label for neuron #i is in the i+1st entry
        p = plot(timevec, data(:, ToPl(iSubPlot) + 2), 'b', ...
            'DisplayName', replace(labels{ToPl(iSubPlot) + 1}, '_', '\_'));
    elseif strcmp(filetype, 'sp')
        % Property #i is in the i+1st column
        p = plot(timevec, data(:, ToPl(iSubPlot) + 1), 'b', ...
            'DisplayName', replace(labels{ToPl(iSubPlot)}, '_', '\_'));
    end
    xlim([xlim1, xlim2]);
    if strcmp(filetype, 'v')
        ylim([-120, 60]);
    elseif strcmp(filetype, 'cli')
    elseif strcmp(filetype, 'sp')
    end

    % Remove X Tick Labels except for the last subplot
    if iSubPlot < nsubplots
        set(gca,'XTickLabel',[])
    end

    % Add stimulation marks
    ax = gca;
    line([stimStartPlot, stimStartPlot], [ax.YLim(1), ax.YLim(2)], ...
        'Color', 'r', 'LineStyle', '--');            % line for stimulation on
    line([stimStartPlot + stimDurPlot, stimStartPlot + stimDurPlot], [ax.YLim(1), ax.YLim(2)], ...
        'Color', 'r', 'LineStyle', '--');            % line for stimulation off

    % For the voltage plot (subplot 1) of special neurons, 
    %   add spike times if available
    if strcmp(filetype, 'sp') && ~isempty(indThis) && iSubPlot == 1
        plot(spiketimesThis, voltageThis, 'r*', 'MarkerSize', 3);
    end

    % Create legend labeling only the trace
    legend([p], 'Location', 'northeast');

    % Add a stimulation line & title for the first subplot and an x-axis label for the last subplot
    if iSubPlot == 1
        if strcmp(filetype, 'v')
            title(['Somatic voltage (mV) for ', replace(filename, '_', '\_')]);
        elseif strcmp(filetype, 'cli')
            title(['Chloride concentration (mM) for ', replace(filename, '_', '\_')]);
        elseif strcmp(filetype, 'sp')
            title(['Traces for ', replace(filename, '_', '\_')]);
        end
    elseif iSubPlot == nsubplots
        xlabel(timelabel);
    end
end

% Save figure
save_all_figtypes(h, figName, figtypes);
%close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [timevec, timelabel, xlim1, xlim2, ...
            stimStartPlot, stimDurPlot, spiketimevec] = ...
                set_time_units(data, tStart, tStop, stimStart, stimDur, spikes)

% Change units of time axis from ms to s if the total time is > 10 seconds
if tStop > 10000
    timevec = data(:, 1)/1000;
    timelabel = 'Time (s)';
    xlim1 = tStart/1000;
    xlim2 = tStop/1000;
    stimStartPlot = stimStart/1000;
    stimDurPlot = stimDur/1000;
    if nargin >= 6
        spiketimevec = spikes(:, 2)/1000;
    else
        spiketimevec = [];
    end
else
    timevec = data(:, 1);
    timelabel = 'Time (ms)';
    xlim1 = tStart;
    xlim2 = tStop;
    stimStartPlot = stimStart;
    stimDurPlot = stimDur;
    if nargin >= 6
        spiketimevec = spikes(:, 2);
    else
        spiketimevec = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [climits, spikes] = get_aux (filetype, filename, inFolder)

% Based on filetype, find range of data values to plot 
%   and spike time data filename
if strcmp(filetype, 'v')
    climits = [-100, 50];
    spifilename = replace(filename, '.singv', '.spi');
elseif strcmp(filetype, 'cli')
    climits = [0, 100];
    spifilename = replace(filename, '.singcli', '.spi');
elseif strcmp(filetype, 'sp')
    climits = [];

    % Start with something like TC[49]_gincr_7.5.singsp
    % First, replace the extension
    temp1 = replace(filename, '.singsp', '.spi');

    % Next, split by ']', which yields a cell array
    temp2 = strsplit(temp1, ']');

    % Save the second part
    tail = temp2{2};

    % Split the first part by '['
    temp3 = strsplit(temp2{1}, '[');

    % Save the first part of that
    head = temp3{1};

    % Finally, concatenate
    spifilename = [head, tail];
end

% Load spike data
spikes = load(fullfile(inFolder, spifilename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

files = dir(fullfile(inFolder, '*.sing*'));
nFiles = numel(files);
fileName = files(i).name;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
