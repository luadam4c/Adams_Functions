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
baseWindowSeconds = 10;
yAmountToStagger = [];
plotIndividualFlag = true;
plotPopulationFlag = false;

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

% Convert to minutes
% TODO: convert_units.m
timeVecsMin = cellfun(@(x) x / S_PER_MIN, timeVecsSec, 'UniformOutput', false);

% Extract the deltaF over F data, ordered from highest to lowest maximum
dFF0ValuesSorted = ...
    cellfun(@(x, y) extract_dFF0_from_slidebook_output(x, y, baseWindowSeconds), ...
            traceTables, num2cell(siSeconds), 'UniformOutput', false);

%% Compute population data
% Compute the average traces
dFF0Means = cellfun(@(x) compute_combined_trace(x, 'mean'), ...
                    dFF0ValuesSorted, 'UniformOutput', false);

% Compute the upper and lower confidence bounds of traces
dFF0Upper95s = cellfun(@(x) compute_combined_trace(x, 'upper95'), ...
                    dFF0ValuesSorted, 'UniformOutput', false);

dFF0Lower95s = cellfun(@(x) compute_combined_trace(x, 'lower95'), ...
                    dFF0ValuesSorted, 'UniformOutput', false);


%% Plot individual traces
if plotIndividualFlag
    figs(1:nFiles) = cellfun(@(x, y, z, w) plot_individual_imaging_traces(x, y, ...
                                            yAmountToStagger, z, w, ...
                                            figTypes, otherArguments), ...
                            timeVecsMin, dFF0ValuesSorted, ...
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
idxTimeVar = find_in_strings(timeStr, varNames, 'SearchMode', 'substrings', ...
                            'MaxNum', 1);

% Extract time vectors in milliseconds
timeVecMs = traceTable{:, idxTimeVar};

% Convert to seconds
timeVecSec = timeVecMs / MS_PER_S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dFF0ValuesSorted = extract_dFF0_from_slidebook_output (traceTable, ...
                                                    siSeconds, baseWindowSeconds)

%% Hard-coded parameters
compStr = 'FlashComp';
meanStr = '_Mean';

%% Remove all variables from the time composite
% Get all the variable names
varNames = traceTable.Properties.VariableNames;

% Find all variable names with the meanStr
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

% Extract the mean values
meanValues = traceTable{:, indMeanVars};

%% Compute the deltaF/F0 for each trace
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

%% Filter out the noise
% Determine the low pass filter cutoff
lowpassCutoff = 1 / (32 * siSeconds);

% Filter out the higher frequencies
dFF0Values = freqfilter(dFF0ValuesRaw, lowpassCutoff, siSeconds, 'FilterType', 'low');

%% Reorder them according to maximum value
% Compute the maximum of values for each vector
maxValues = compute_stats(dFF0Values, 'max');

% Sort the values from high to low
[~, origIndex] = sort(maxValues, 'descend');

% TODO: Reorder cell numbers

% Reorder the traces
dFF0ValuesSorted = dFF0Values(:, origIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig = plot_individual_imaging_traces (timeVecMin, dFF0Values, ...
                                            yAmountToStagger, figTitle, ...
                                            figName, figTypes, otherArguments)

%% Plot the mean values
% Decide on the amount to stagger on the y axis
if isempty(yAmountToStagger)
    % Compute the range of values for each vector
    rangeValues = compute_stats(dFF0Values, 'range');

    % Use 1/4 of the average range
    yAmountToStagger = mean(rangeValues) / 4;
end

% Create figure
fig = set_figure_properties('AlwaysNew', true, 'Units', 'inches', ...
                            'Width', 1.5, 'Height', 2.3);

xLabel = 'Time (min)';
yLabel = 'Cell Number';
                        
% Plot as a staggered plot
[fig, ax] = ...
    plot_traces(timeVecMin, dFF0Values, 'PlotMode', 'staggered', ...
                'YAmountToStagger', yAmountToStagger, 'YBase', 0, ...
                'XLabel', xLabel, 'YLabel', yLabel, ...
                'LineWidth', 1, 'ColorMap', 'k', ...
                'FigTitle', figTitle, 'FigHandle', fig, otherArguments);

% Update the figure to look pretty
fig = update_figure_for_corel(fig);

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