function figs = plot_calcium_imaging_traces (varargin)
%% Plots calcium imaging traces
% Usage: figs = plot_calcium_imaging_traces (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       figs        - figures plotted
%                   specified as a Figure object handle array
%
% Arguments:
%       varargin    - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == {'png', 'epsc'}
%                   - Any other parameter-value pair for 
%                       plot_traces() or plot_tuning_curve()
%
% Requires:
%       cd/all_files.m
%       cd/compute_stats.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/find_first_match.m
%       cd/find_in_strings.m
%       cd/freqfilter.m
%       cd/isfigtype.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%       cd/save_all_figtypes.m
%       cd/struct2arglist.m
%       cd/update_figure_for_corel.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-19 Created by Adam Lu
% 2019-09-30 Made into a function
% 

%% Hard-coded parameters
S_PER_MIN = 60;

% TODO: Make optional arguments
baseWindowSeconds = 1; %10;
yAmountToStagger = [];
plotIndividualFlag = true;
plotPopulationFlag = false;
toNormalizeByBaseline = false; %true;
toFilter = false; %true;
toSort = false; %true;

%% Default values for optional arguments
figTypesDefault = {'png', 'epsc'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for plot_traces() or plot_tuning_curve()
otherArguments = iP.Unmatched;

%% Preparation
% Get all paths to trace text files
[~, textPaths] = all_files('Suffix', 'time_series', 'Extension', 'txt', ...
                            'ForceCellOutput', true);

% Count the number of files
nFiles = numel(textPaths);

% Initialize figure array
figs = gobjects(nFiles, 1);

% Get the base names
textBases = extract_fileparts(textPaths, 'base');

% Create figure titles for individual plots
figTitlesSingles = strcat('All traces for ', replace(textBases, '_', '\_'));

% Create figure names for individual plots
figNamesSingles = strcat(textBases, '_individual.png');

% Extract the common prefix
commonPrefix = extract_fileparts(textPaths, 'commonprefix');

% Create a figure title for the population plot
figTitlePop = strcat('Averaged traces for', replace(commonPrefix, '_', '\_'));

% Create a figure name for the population plot
figNamePop = strcat(commonPrefix, '_population.png');

%% Read in and process data from each file
% Read in the table(s)
traceTables = cellfun(@readtable, textPaths, 'UniformOutput', false);

% Extract time vector(s) in seconds
timeVecsSec = cellfun(@extract_time_vector_from_slidebook_output, ...
                        traceTables, 'UniformOutput', false);

% Compute the sampling interval(s) in seconds
siSeconds = compute_sampling_interval(timeVecsSec, 'IsRegular', false);

% Extract the deltaF over F data, ordered from highest to lowest maximum
[dFF0ValuesSorted, traceLabels, timeVecsSec] = ...
    cellfun(@(x, y, z) extract_dFF0_from_slidebook_output(x, y, z, ...
                                    baseWindowSeconds, toNormalizeByBaseline, ...
                                    toFilter, toSort), ...
            traceTables, timeVecsSec, num2cell(siSeconds), 'UniformOutput', false);

% Convert to minutes
% TODO: convert_units.m
timeVecsMin = cellfun(@(x) x / S_PER_MIN, timeVecsSec, 'UniformOutput', false);

%% Compute population data
if plotPopulationFlag
    % Compute the average traces
    dFF0Means = cellfun(@(x) compute_combined_trace(x, 'mean'), ...
                        dFF0ValuesSorted, 'UniformOutput', false);

    % Compute the upper and lower confidence bounds of traces
    dFF0Upper95s = cellfun(@(x) compute_combined_trace(x, 'upper95'), ...
                        dFF0ValuesSorted, 'UniformOutput', false);

    dFF0Lower95s = cellfun(@(x) compute_combined_trace(x, 'lower95'), ...
                        dFF0ValuesSorted, 'UniformOutput', false);
end

%% Plot individual traces
if plotIndividualFlag
    figs(1:nFiles) = cellfun(@(x, y, z, w, v) ...
                                plot_individual_imaging_traces(x, y, z, ...
                                            yAmountToStagger, w, v, ...
                                            figTypes, toNormalizeByBaseline, otherArguments), ...
                            timeVecsMin, dFF0ValuesSorted, traceLabels, ...
                            figTitlesSingles, figNamesSingles);
end

%% Plot population data
if plotPopulationFlag
    figs(nFiles + 1) = plot_averaged_imaging_traces(timeVecsMin, dFF0Means, ...
                            dFF0Upper95s, dFF0Lower95s, ...
                            figTitlePop, figNamePop, figTypes, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timeVecSec = extract_time_vector_from_slidebook_output (traceTable)

%% Hard-coded parameters
MS_PER_S = 1000;
timeStr = '_ms';

% Get all the variable names
varNames = traceTable.Properties.VariableNames;

% Find all variable names with the timeStr
idxTimeVar = find_first_match(timeStr, varNames, ...
                            'MatchMode', 'parts', 'IgnoreCase', true);

% Extract time vectors in milliseconds
timeVecMs = traceTable{:, idxTimeVar};

% Convert to seconds
timeVecSec = timeVecMs / MS_PER_S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dFF0Values, traceLabels, timeVecsSec] = ...
                extract_dFF0_from_slidebook_output (traceTable, timeVecsSec, ...
                                    siSeconds, baseWindowSeconds, ...
                                    toNormalizeByBaseline, toFilter, toSort)

%% Hard-coded parameters
compStr = 'FlashComp';
meanStr = '_Mean';

% TODO: Make this optional
channelsToCombine = {'488_561Dual', '561Prime'}; % {};

%% Remove all variables from the time composite
% Get all the variable names
varNames = traceTable.Properties.VariableNames;

% Find all variable names with the compStr
[~, compVarNames] = find_in_strings(compStr, varNames, ...
                                'SearchMode', 'substrings');

% Remove the variables
if ~isempty(compVarNames)
    traceTable = removevars(traceTable, compVarNames);
end

%% Extract the mean F values
% Get all the variable names
varNames = traceTable.Properties.VariableNames;

% Find all variable names with the meanStr
[indMeanVars, meanVarNames] = find_in_strings(meanStr, varNames, ...
                                'SearchMode', 'substrings');

% TODO: extract the cell numbers from meanVarNames
% Create default trace labels
traceLabels = meanVarNames;

% Extract the mean values
meanValues = traceTable{:, indMeanVars};

%% Extract endpoints if requested
% TODO: Make optional
% timeWindowSeconds = [];
timeWindowSeconds = timeVecsSec(1) + [20, 70];

if ~isempty(timeWindowSeconds)
    endPoints = find_window_endpoints(timeWindowSeconds, timeVecsSec);

    meanValues = extract_subvectors(meanValues, 'EndPoints', endPoints);
    timeVecsSec = extract_subvectors(timeVecsSec, 'EndPoints', endPoints);
end

%% Combine mean values if requested
% TODO: Make this optional
% Count the number of channels to combine
nChannelsToCombine = numel(channelsToCombine);

% 
if nChannelsToCombine > 0 && nChannelsToCombine <= 2
    % Get the name of the first channel
    channelFirst = channelsToCombine{1};

    % Find all variables for this channel
    [indFirstVars, namesFirstVars] = ...
        find_in_strings(channelFirst, meanVarNames, 'SearchMode', 'substrings');
                                
    % Extract trace labels
    traceSuffixes = extractAfter(namesFirstVars, channelFirst);

    % Get the name of the second channel
    channelNext = channelsToCombine{2};

    % Find all variables for this channel
    [indSecondVars, namesSecondVars] = ...
        cellfun(@(x) find_in_strings({channelNext, x}, meanVarNames, ...
                    'SearchMode', 'substrings'), traceSuffixes, 'UniformOutput', false);
    indSecondVars = cell2num(indSecondVars);

    % Combine the traces in turn
    combinedTraces = arrayfun(@(x, y) sum(meanValues(:, [x, y]), 2), ...
            indFirstVars, indSecondVars, 'UniformOutput', false);

    % Force as a matrix
    meanValues = force_matrix(combinedTraces);
    
    % TODO: Make this better
    traceLabels = extractAfter(extractBefore(traceSuffixes, meanStr), '_');
elseif nChannelsToCombine > 0
    error('Not implemented yet!');
end

%% Compute the deltaF/F0 for each trace
if toNormalizeByBaseline
    % Convert the baseline window to samples
    baseWindowSamples = [1, ceil(baseWindowSeconds / siSeconds)];

    % Extract baseline values
    baseValues = extract_subvectors(meanValues, 'EndPoints', baseWindowSamples);

    % Compute the initial values F0
    initValues = compute_stats(baseValues, 'mean', 'IgnoreNan', true);

    % Match the number of rows
    initValues = match_row_count(initValues, size(meanValues, 1));

    % Compute the deltaF/F0 = (F - F0)/F0
    dFF0ValuesRaw = (meanValues - initValues) ./ initValues;
else
    % Use F instead
    dFF0ValuesRaw = meanValues;
end

%% Filter out the noise
if toFilter
    % Determine the low pass filter cutoff
    lowpassCutoff = 1 / (32 * siSeconds);

    % Filter out the higher frequencies
    dFF0Values = freqfilter(dFF0ValuesRaw, lowpassCutoff, siSeconds, 'FilterType', 'low');
else
    dFF0Values = dFF0ValuesRaw;
end

%% Reorder them according to maximum value
if toSort
    % Compute the maximum of values for each vector
    maxValues = compute_stats(dFF0Values, 'max');

    % Sort the values from high to low
    [~, origIndex] = sort(maxValues, 'descend');

    % TODO: Reorder cell numbers

    % Reorder the traces
    dFF0Values = dFF0Values(:, origIndex);

    % Reorder trace labels
    traceLabels = traceLabels(origIndex);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig = plot_individual_imaging_traces (timeVecMin, dFF0Values, ...
                                    traceLabels, yAmountToStagger, figTitle, ...
                                    figName, figTypes, toNormalizeByBaseline, otherArguments)

% Decide on the plotting mode
if toNormalizeByBaseline
    plotMode = 'staggered';
    yLabel = 'Cell Number';
else
    plotMode = 'overlapped';
    yLabel = 'Absolute F value';
end

%% Plot the mean values
% Decide on the amount to stagger on the y axis
if isempty(yAmountToStagger)
    % Compute the range of values for each vector
    rangeValues = compute_stats(dFF0Values, 'range');

    % Use 1/4 of the average range
    yAmountToStagger = mean(rangeValues) / 4;
end

% Create figure
fig = set_figure_properties('AlwaysNew', true);

xLabel = 'Time (min)';
                        
% Plot as a staggered plot
plot_traces(timeVecMin, dFF0Values, 'PlotMode', plotMode, ...
                'YAmountToStagger', yAmountToStagger, 'YBase', 0, ...
                'XLabel', xLabel, 'YLabel', yLabel, ...
                'TraceLabels', traceLabels, ...
                'LineWidth', 1, ...
                'LegendLocation', 'eastoutside', ...
                'FigTitle', figTitle, 'FigHandle', fig, otherArguments);
% 'ColorMap', 'k', 

% Save figure
save_all_figtypes(fig, figName, figTypes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig = plot_averaged_imaging_traces (timeVecsMin, dFF0Means, ...
                    dFF0Upper95s, dFF0Lower95s, figTitlePop, figNamePop, ...
                    figTypes, otherArguments)

%% Preparation
% Count the number of slices
nSlices = numel(dFF0Means);

% Create figure
fig = set_figure_properties('AlwaysNew', true, 'Units', 'inches', ...
                            'Width', 1.5, 'Height', 1.5);

% Create subplots
[fig, ax] = create_subplots(nSlices, 1, 'FigHandle', fig);

% Colors
% TODO: Make these optional arguments
colorMap = decide_on_colormap([], nSlices);

% Plot as subplots
for iSlice = 1:nSlices
    timeVec = timeVecsMin{iSlice};
    mean = dFF0Means{iSlice};
    lower95 = dFF0Lower95s{iSlice};
    upper95 = dFF0Upper95s{iSlice};

    if iSlice == 1
        yLimits = [-0.2, 0.4];
    elseif iSlice == 2
        yLimits = [-1, 1.4];
    elseif iSlice == 3
        yLimits = [-0.06, 0.12];
    end
       
    axes(ax(iSlice));
    plot_tuning_curve(timeVec, mean, 'LowerCI', lower95, ...
                        'UpperCI', upper95, 'ReadoutLimits', yLimits, ...
                        'FigTitle', 'suppress', ...
                        'PLabel', 'suppress', 'ReadoutLabel', 'suppress', ...
                        'ColorMap', colorMap(iSlice, :), otherArguments);
    
end

% Update the figure to look pretty
fig = update_figure_for_corel(fig);

% Save figure
save_all_figtypes(fig, figNamePop, figTypes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
