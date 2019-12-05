% m3ha_plot_figure02.m
%% Plots Figure 02 for the GAT Blocker paper
%
% Requires:
%       cd/all_files.m
%       cd/apply_iteratively.m
%       cd/argfun.m
%       cd/copy_into.m
%       cd/decide_on_colormap.m
%       cd/extract_columns.m
%       cd/extract_fileparts.m
%       cd/extract_substrings.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_create_cell_info_table.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_load_sweep_info.m
%       cd/plot_scale_bar.m
%       cd/plot_traces.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m

% File History:
% 2019-11-25 Created by Adam Lu

%% Hard-coded parameters
% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
matFilesDir = fullfile(parentDirectory, 'data_dclamp', 'take4', 'matfiles');

% Files
datalogPath = fullfile(figure02Dir, 'dclampdatalog_take4.csv');
cellNamePattern = '[A-Z][0-9]{6}';

% Flags
saveCellInfo = false;
copyExampleFiles = false; %true;
plotExamplesFlag = true;
plotViolinPlotsFlag = false; %true;
plotBarPlotsFlag = false; %true;

% Analysis settings
% Note: must be consistent with m3ha_compute_statistics.m
measuresOfInterest = {'ltsAmplitude'; 'ltsMaxSlope'; ...
                    'ltsConcavity'; 'ltsProminence'; ...
                    'ltsWidth'; 'ltsOnsetTime'; 'ltsTimeJitter'; ...
                    'ltsProbability'; 'spikesPerLts'; ...
                    'spikeMaxAmp'; 'spikeMinAmp'; ...
                    'spikeFrequency'; 'spikeAdaptation'
                    'burstOnsetTime'; 'burstTimeJitter'; ...
                    'burstProbability'; 'spikesPerBurst'};

% Plot settings
exampleGincr = 200;
exampleHeight = 9;             % in centimeters
exampleWidth = 8.5;             % in centimeters
exampleLineWidth = 0.5;
exampleXlimits = [800, 2200];
exampleYlimits = [-100, 10];
exampleYTicks = [-80, -50, -20];

pharmAll = [1; 2; 3; 4];          
pharmLabels = {'{\it d}-Control', '{\it d}-GAT1 Block', ...
                '{\it d}-GAT3 Block', '{\it d}-Dual Block'};
pharmLabelsShort = {'{\it d}-Con', '{\it d}-GAT1', ...
                    '{\it d}-GAT3', '{\it d}-Dual'};
gIncrAll = [25; 50; 100; 200; 400; 800];
gIncrLabels = {'25%', '50%', '100%', '200%', '400%', '800%'};
conditionLabel2D = 'pharm_1-4_gincr_200';
pCond2D = num2cell(pharmAll);
gCond2D = 200;
conditionLabel3D = 'pharm_1-4_gincr_all';
pCond3D = num2cell(pharmAll);
gCond3D = num2cell(gIncrAll);

violinFigHeight = 5;            % in centimeters
violinFigWidth = 3.4;           % in centimeters
violinRelativeBandWidth = 0.1;  % bandwidth relative to data range
medianColor = rgb('GreenYellow');     % color of median circle
medianSize = 6;                % size of median circle in points
bar3FigHeight = 6;              % in centimeters
bar3FigWidth = 6;               % in centimeters

figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load sweep info
% Read from datalogPath
swpInfo = m3ha_load_sweep_info('Directory', figure02Dir);

% Look for all import logs
[~, exampleLogPaths] = all_files('Directory', figure02Dir, ...
                        'Suffix', 'imported_files', 'Extension', 'txt');

% Extract example cell names
exampleCellNames = ...
    extract_substrings(exampleLogPaths, 'Regexp', cellNamePattern);

%% Create cell info table
if saveCellInfo
    cellInfo = m3ha_create_cell_info_table('SwpInfo', swpInfo);
end

%% Copy example traces
if copyExampleFiles
    % Copy files from the log files
    cellfun(@(x) copy_files_from_log(x, matFilesDir, figure02Dir), ...
            exampleLogPaths);
end

%% Plot example traces
if plotExamplesFlag
    % Plot example traces for each cell
    handles = cellfun(@(x) m3ha_plot_example_traces(x, figure02Dir, ...
                            swpInfo, exampleGincr, ...
                            exampleXlimits, exampleYlimits, ...
                            exampleYTicks, pharmLabels, exampleLineWidth, ...
                            exampleHeight, exampleWidth, figTypes), ...
                        exampleCellNames);

end

%% Plot 2D violin plots
if plotViolinPlotsFlag
    % Construct stats table path
    statsPath2D = fullfile(figure02Dir, strcat(conditionLabel2D, '_stats.mat'));

    % Load or compute statistics
    if isfile(statsPath2D)
        % Load stats table
        disp('Loading statistics for 2D violin plots ...');
        load(statsPath2D, 'statsTable');
    else
        % Compute statistics for all features
        disp('Computing statistics for 2D violin plots ...');
        statsTable = m3ha_compute_statistics('PharmConditions', pCond2D, ...
                                                'GIncrConditions', gCond2D);

        % Save stats table
        save(statsPath2D, 'statsTable', '-v7.3');
    end

    % Restrict to measures of interest
    statsTable = statsTable(measuresOfInterest, :);

    % Extract variables
    allMeasureTitles = statsTable.measureTitle;
    allMeasureStrs = statsTable.measureStr;
    allValues = statsTable.allValues;

    % Create figure bases
    allFigBases2D = combine_strings({allMeasureStrs, conditionLabel2D});

    % Create full path bases
    allFigPathBases2D = fullfile(figure02Dir, allFigBases2D);

    % Plot all 2D violin plots
    disp('Plotting 2D violin plots ...');
    % cellfun(@(a, b, c) m3ha_plot_violin(a, violinRelativeBandWidth, b, ...
    %                             pharmLabels, c, ...
    %                             violinFigHeight, violinFigWidth, figTypes), ...
    %         allValues, allMeasureTitles, allFigPathBases2D);
    cellfun(@(a, b, c) m3ha_plot_violin(a, violinRelativeBandWidth, ...
                                medianColor, medianSize, b, ...
                                pharmLabelsShort, c, ...
                                violinFigHeight, violinFigWidth, figTypes), ...
            allValues, allMeasureTitles, allFigPathBases2D);
end

%% Plot 3D bar plots
if plotBarPlotsFlag
    % Construct stats table path
    statsPath3D = fullfile(figure02Dir, strcat(conditionLabel3D, '_stats.mat'));

    % Load or compute statistics
    if isfile(statsPath3D)
        % Load stats table
        disp('Loading statistics for 3D bar plots ...');
        load(statsPath3D, 'statsTable');
    else
        % Compute statistics for all features
        disp('Computing statistics for 3D bar plots ...');
        statsTable = m3ha_compute_statistics('PharmConditions', pCond3D, ...
                                                'GIncrConditions', gCond3D);

        % Save stats table
        save(statsPath3D, 'statsTable', '-v7.3');
    end

    % Restrict to measures of interest
    statsTable = statsTable(measuresOfInterest, :);

    % Extract variables
    allMeasureTitles = statsTable.measureTitle;
    allMeasureStrs = statsTable.measureStr;
    allMeanValues = statsTable.meanValue;
    allUpper95Values = statsTable.upper95Value;

    % Create figure bases
    allFigBases3D = combine_strings({allMeasureStrs, conditionLabel3D});

    % Create full path bases
    allFigPathBases3D = fullfile(figure02Dir, allFigBases3D);

    % Plot all 3D bar plots
    disp('Plotting 3D bar plots ...');
    cellfun(@(a, b, c, d) m3ha_plot_bar3(a, b, c, ...
                    pharmLabels, gIncrLabels, ...
                    d, bar3FigHeight, bar3FigWidth, figTypes), ...
            allMeanValues, allUpper95Values, ...
            allMeasureTitles, allFigPathBases3D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_files_from_log(logFile, dirFrom, dirTo)

% Read the log file
logTable = readtable(logFile, 'ReadVariableNames', false, 'Delimiter', ' ');

% Extract .mat file names
matFiles = logTable.Var1;

% Construct full paths to .mat files
matPaths = fullfile(dirFrom, matFiles);

% Copy files
copy_into(matPaths, dirTo);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_example_traces (cellName, figure02Dir, ...
                                            swpInfo, gIncrOfInterest, ...
                                            xLimits, yLimits, yTicks, ...
                                            pharmLabels, lineWidth, ...
                                            figHeight, figWidth, figTypes)
%% Plots example traces
% TODO: Pull out as its own function with only cellName as the required argument

NS_PER_US = 1000;

%% Hard-coded parameters
conductanceLabel = 'Conductance (nS)';
conductanceYLimits = [0, 10];
timeBarLength = 100;
timeBarUnits = 'ms';

% Look for all .mat files
[~, matPathsAll] = all_files('Directory', figure02Dir, ...
                        'KeyWord', cellName, 'Extension', 'mat');

% Extract file bases
matBasesAll = extract_fileparts(matPathsAll, 'base');

% Get all the conductance amplitudes
gIncrAll = swpInfo{matBasesAll, 'grow'};

% Restrict to the conductance amplitude of interest
exampleBases = matBasesAll(gIncrAll == gIncrOfInterest);

% Get all the pharm conditions
pharmAll = swpInfo{exampleBases, 'prow'};

% Sort according to pharm condition
[~, origInd] = sort(pharmAll);
exampleBases = exampleBases(origInd);

% Import traces
data = m3ha_import_raw_traces(exampleBases, 'SweepInfoAll', swpInfo, ...
                                'OutFolder', figure02Dir);

% Extract time, voltage and conductance vectors
[tVecs, vVecs, gVecs] = extract_columns(data, [1:2, 4]);

% Create figure for voltage traces
figV = set_figure_properties('AlwaysNew', true);

% Plot the voltage traces
tracesV = plot_traces(tVecs, vVecs, ...
            'PlotMode', 'parallel', 'LineWidth', lineWidth, ...
            'LinkAxesOption', 'xy', 'FigTitle', 'suppress', ...
            'XLimits', xLimits, 'Ylimits', yLimits, ...
            'XLabel', 'suppress', 'YLabel', pharmLabels, ...
            'LegendLocation', 'suppress', ...
            'FigHandle', figV);

% Get voltage subplots
axV = tracesV.subPlots;

% Remove x axis
set(axV, 'XTick', []);

% Plot a time bar
plot_scale_bar('x', 'BarLength', timeBarLength, 'BarUnits', timeBarUnits);

% Update figure for CorelDraw
update_figure_for_corel(figV, 'Units', 'centimeters', ...
                        'YTickLocs', yTicks, ...
                        'Height', figHeight, 'Width', figWidth);

% Create figure base
figBase = sprintf('%s_gincr_%g_examples', cellName, gIncrOfInterest);
figPathBase = fullfile(figure02Dir, figBase);

% Save the figure
save_all_figtypes(figV, figPathBase, figTypes);

% Create figure for conductance traces
%   Note: prevent the legend from updating
figG = set_figure_properties('AlwaysNew', true, 'defaultLegendAutoUpdate', 'off');

% Convert conductance to nS
% TODO: use convert_units.m
gVecs = cellfun(@(x) x .* NS_PER_US, gVecs, 'UniformOutput', false);

% Plot the conductance traces
tracesG = plot_traces(tVecs, gVecs, ...
            'PlotMode', 'overlapped', 'LineWidth', lineWidth, ...
            'FigTitle', 'suppress', ...
            'XLimits', xLimits, 'YLimits', conductanceYLimits, ...
            'XLabel', 'suppress', 'YLabel', conductanceLabel, ...
            'TraceLabels', pharmLabels, ...
            'LegendLocation', 'northeast', ...
            'FigHandle', figG);

% Get conductance trace axes
axG = tracesG.subPlots;

% Remove x axis
set(axG, 'XTick', []);

% Plot a time bar
plot_scale_bar('x', 'BarLength', timeBarLength, 'BarUnits', timeBarUnits);

% Update figure for CorelDraw
update_figure_for_corel(figG, 'Units', 'centimeters', ...
                        'Height', figHeight / 4, 'Width', figWidth);

% TODO: match_axes_size.m 
%   match_axes_size(axG, axV(1), 'width')

% Get the axes width for the voltage plot
voltageSubPlotPosition = get(axV(1), 'Position');

% Update axes width to be consistent with the voltage plot
set_axes_properties('AxesHandle', axG, 'Width', voltageSubPlotPosition(3));

% Create figure base
figBase = sprintf('%s_gincr_%g_conductance', cellName, gIncrOfInterest);
figPathBase = fullfile(figure02Dir, figBase);

% Save the figure
save_all_figtypes(figG, figPathBase, figTypes);

% Output figure handles
handles.figG = figG;
handles.figV = figV;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3ha_plot_violin (allValues, violinRelativeBandWidth, ...
                            medianColor, medianSize, ...
                            measureTitle, pharmLabels, ...
                            figPathBase, figHeight, figWidth, figTypes)

% Hard-coded parameters
MS_PER_S = 1000;
xTickAngle = 320;

% Create figure for conductance traces
fig = set_figure_properties('AlwaysNew', true);

% Count the number of groups
nGroups = numel(pharmLabels);

% Convert onset times from ms to seconds
%{
if contains(measureTitle, 'onset')
    % Update values
    allValues = cellfun(@(x) x ./ MS_PER_S, allValues, 'UniformOutput', false);

    % Update title
    measureTitle = replace(measureTitle, 'ms', 's');
end
%}

% Decide on the color map
cm = decide_on_colormap([], nGroups);

% Set the color map
%   TODO: Apply this in violinplot?
colormap(cm);

% TODO: plot_violin.m
% TODO: plot_jitter.m

% Force as a numeric array
allValues = force_matrix(allValues);

% Compute range of all values
rangeValues = apply_iteratively(@max, allValues) - apply_iteratively(@min, allValues);

% Compute the bandwidth for the kernel density estimates
bandWidth = violinRelativeBandWidth * rangeValues;

% Plot a violin plot
% violinplot(allValues, pharmLabels);
violinplot(allValues, pharmLabels, 'BandWidth', bandWidth);

% Plot the data points for each cell
% plotSpread(allValues);
% Set x tick labels
% xticklabels(pharmLabels);

% Modify x limits
xlim([0.5, nGroups + 0.5]);

% Modify x tick angle
xtickangle(xTickAngle);

% Set y label
ylabel(measureTitle);

% Find all median scatters and make the face color green and size bigger
medianScatters = findobj(gca, 'Type', 'Scatter', ...
                        'MarkerEdgeColor', [0.5, 0.5, 0.5]);
set(medianScatters, 'MarkerFaceColor', medianColor);

% Save the figure
save_all_figtypes(fig, [figPathBase, '_orig'], 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth, ...
                        'ScatterMarkerSize', 3);

% Fix axes position
set(gca, 'Position', [0.2356, 0.1947, 0.6694, 0.7303]);

% Update median size
set(medianScatters, 'SizeData', medianSize^2);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3ha_plot_bar3(meanValue, upper95Value, ...
                        measureTitle, pharmLabels, gIncrLabels, ...
                        figPathBase, figHeight, figWidth, figTypes)

% Create figure for conductance traces
fig = set_figure_properties('AlwaysNew', true);

% Flip the g incr axis
[meanValue, upper95Value, gIncrLabels] = ...
    argfun(@fliplr, meanValue, upper95Value, gIncrLabels);

% Set x and y tick labels
xTickLabels = pharmLabels;
yTickLabels = gIncrLabels;

% TODO: Add the following to plot_bar.m?
% Hard-coded parameters
relativeBarWidth = 0.2;
xTickAngle = 320;
barSeparation = 1;

% Decide on the color map
cm = decide_on_colormap([], 4);

% Set the color map
colormap(cm);

% Prepare for bar3
meanValueTransposed = transpose(meanValue);
upper95ValueTransposed = transpose(upper95Value);

% Plot the means as bars
bar3(meanValueTransposed, relativeBarWidth, 'detached');

% Plot error bars
% TODO: Incorporate into plot_error_bar.m?

% Set the relative error bar width to be the same as the bars themselves
%   Note: error bar width must not exceed the bar width, 
%           otherwise the edges would be cut off
relativeErrorBarWidth = relativeBarWidth;

% Compute the actual error bar width
errorBarWidth = relativeErrorBarWidth * barSeparation;

% Compute the x and y values corresponding to each data point
[xValues, yValues] = meshgrid(1:numel(xTickLabels), 1:numel(yTickLabels));

% Compute the left and right positions of the horizontal parts of the error bars
xPosBarLeft = xValues - errorBarWidth / 2;
xPosBarRight = xValues + errorBarWidth / 2;

% Plot the vertical part of the error bars
arrayfun(@(a, b, c, d, e, f) line([a, b], [c, d], [e, f], 'Color', 'k'), ...
            xValues, xValues, yValues, yValues, ...
            meanValueTransposed, upper95ValueTransposed);

% Plot the horizontal part of the error bars
arrayfun(@(a, b, c, d, e, f) line([a, b], [c, d], [e, f], 'Color', 'k'), ...
            xPosBarLeft, xPosBarRight, yValues, yValues, ...
            upper95ValueTransposed, upper95ValueTransposed);

% Plot z axis label
zlabel(measureTitle);

% Set x tick labels
set(gca, 'XTickLabel', xTickLabels);

% Set x tick angle
xtickangle(xTickAngle);

% Set y tick labels
set(gca, 'YTickLabel', yTickLabels);

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

barColor = 'k';
barLineWidth = 1;
xPosNormalized = 0.8;
yPosBarNormalized = 0.1;
yPosTextNormalized = 0.05;
xPosBar = xLimits(1) + xPosNormalized * diff(xLimits);
yPosBar = yLimits(1) + yPosBarNormalized * diff(yLimits);
xPosText = xLimits(1) + xPosNormalized * diff(xLimits);
yPosText = yLimits(1) + yPosTextNormalized * diff(yLimits);
xLimitsBar = xPosBar + [0, barLengthMs];
plot_horizontal_line(yPosBar, 'XLimits', xLimitsBar, ...
                    'Color', barColor, 'LineWidth', barLineWidth, ...
                    'AxesHandle', subPlots(end));
text(xPosText, yPosText, barText);

measuresOfInterest = {'ltsOnsetTime'; 'ltsProbability'; 'spikesPerLts'; ...
                    'burstOnsetTime'; 'burstProbability'; 'spikesPerBurst'};

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
