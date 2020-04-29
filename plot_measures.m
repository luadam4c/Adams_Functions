function plot_measures (varargin)
%% Plots all measures of interest across slices
% Usage: plot_measures (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       plot_measures('PlotAll', true);
%       plot_measures('PhaseNumbers', 1:2);
%       plot_measures('SweepsRelToPhase2', -19:30);
%       plot_measures('SweepNumbers', 1:50);
%       plot_measures('PlotChevron', true, 'SweepsRelToPhase2', -9:31, 'ColorMap', 'k');
%
% Arguments:
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'tuning'    - circles
%                       'bar'       - horizontal bars
%                   default == 'tuning'
%                   - 'ComputeChevronFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeTimeTablesFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputePopAverageFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeNormTablesFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeNormChevronFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeNormTimeTablesFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeNormPopAverageFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeSmoothTablesFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeSmoothTimeTablesFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == fPlotType',tuning,'alse
%                   - 'ComputeSmoothNormTimeTablesFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ComputeSmoothNormPopAverageFlag': whether to compute TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotAllFlag': whether to plot everything
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotByFileFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotByPhaseFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotPopAverageFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotChevronFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotNormByFileFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotNormByPhaseFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotNormPopAverageFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotNormChevronFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSmoothByFileFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSmoothNormByFileFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSmoothNormPopAvgFlag': whether to plot TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveOutliersInPlot': whether to remove outliers 
%                                               in time series plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': working directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'InFolder': directory to read files from
%                   must be a string scalar or a character vector
%                   default == same as directory
%                   - 'OutFolder': directory to place output spreadsheet files
%                   must be a string scalar or a character vector
%                   default == fullfile(inFolder, [create_time_stamp, '_', ...
%                                           outDirSuffix])
%                   - 'FigFolder': directory to place figures
%                   must be a string scalar or a character vector
%                   default == same as outFolder
%                   - 'PhaseNumbers': phase numbers to restrict to
%                   must be a numeric vector
%                   default == []
%                   - 'SweepNumbers': sweep numbers to restrict to
%                   must be a numeric vector
%                   default == []
%                   - 'SweepsRelToPhase2': sweep numbers relative to phase 2
%                                           to restrict to
%                   must be a numeric vector
%                   default == []
%                   - 'NSweepsLastOfPhase': number of sweeps at 
%                                           the last of a phase
%                   must be a positive integer scalar
%                   default == 10
%                   - 'NSweepsToAverage': number of sweeps to average
%                   must be a positive integer scalar
%                   default == 5
%                   - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'maxRange2Mean'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be a nonnegative scalar
%                   default == 40%
%                   - 'SweepLengthSec': sweep length in seconds
%                   must be a nonnegative scalar
%                   default == 60 seconds
%                   - 'TimeLabel': time axis label
%                   must be a string scalar or a character vector
%                   default == 'Time'
%                   - 'PhaseLabel': phase axis label
%                   must be a string scalar or a character vector
%                   default == 'Phase'
%                   - 'PhaseStrings': phase strings
%                   must be a string vector or a cell array of character vectors
%                   default == {'Baseline', 'Wash-on', 'Wash-out'}
%                   - 'VarsToPlot': variables to plot
%                   must be a string vector or a cell array of character vectors
%                   default == varsToPlotAll
%                   - 'VarLabels': variable labels
%                   must be a string vector or a cell array of character vectors1
%                   default == varLabelsAll
%                   - Any other parameter-value pair 
%                       for plot_table() and plot_tuning_curve()
%
% Requires:
%       cd/argfun.m
%       cd/check_dir.m
%       cd/combine_variables_across_tables.m
%       cd/compute_phase_average.m
%       cd/compute_population_average.m
%       cd/copy_into.m
%       cd/count_vectors.m
%       cd/create_label_from_sequence.m
%       cd/create_indices.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/force_matrix.m
%       cd/is_var_in_table.m
%       cd/ismatch.m
%       cd/modify_table.m
%       cd/plot_chevron.m
%       cd/plot_chevron_bar_inset.m
%       cd/plot_table.m
%       cd/plot_tuning_curve.m
%       cd/renamevars.m
%       cd/save_all_figtypes.m
%       cd/set_default_flag.m
%       cd/set_figure_properties.m
%       cd/unique_custom.m
%       cd/update_figure_for_corel.m
%       cd/write_timetable.m
%
% Used by:
%       cd/clc2_analyze.m
%       cd/Glucose_analyze.m

% File History:
% 2019-03-15 Created by Adam Lu
% 2019-03-25 Now colors by phase number
% 2019-04-08 Renamed as plot_measures.m
% 2019-06-11 Now plots both normalized and original versions of the Chevron Plot
% 2019-08-06 Now plots markers for the Chevron plots
% 2019-08-06 Added phaseNumbers
% 2019-08-06 Renamed '_averaged' -> '_chevron'
% 2019-08-07 Added input parser and plotAllFlag
% 2019-08-07 Added directory, inFolder, outFolder
% 2019-08-07 Extracted specific usage to clc2_analyze.m
% 2019-08-07 Added 'SweepNumbers' as an optional argument
%               to allow the restriction to certain sweep numbers
% 2019-08-09 Now runs both T test and Wilcoxon signed-rank test for p-values
% 2019-08-20 Added 'SweepsRelToPhase2' as an optional argument
% 2019-08-20 Fixed population average computation to ignore NaNs
% 2019-08-20 Added 'RemoveOutliersInPlot' as an optional argument
% 2019-08-20 Now sets 'ReturnLastTrial' to be true for compute_phase_average 
%               to avoid NaNs
% 2019-08-29 Added computation of smoothed tables
% 2020-02-19 Now passes extra arguments to plot_table and plot_tuning_curve
% 2020-02-19 Now plots Chevron plots for CorelDraw by default
% TODO: Plot all plots for CorelDraw by default
% TODO: Restrict to certain time windows 'TimeWindow'
% TODO: Pull out many code to functions

%% Hard-coded parameters
validSelectionMethods = {'notNaN', 'maxRange2Mean'};
validPlotTypes = {'tuning', 'bar'};
outDirSuffix = 'population_measures';
smoothFunc = @(x) movingaveragefilter(x, 5, 1);

%% TODO: Make optional argument
backupSheets = true;
plotForCorel = true;
figTypes = {'png', 'epsc'};

% Must be consistent with parse_multiunit.m
varsToPlotAll = {'oscIndex1'; 'oscIndex2'; 'oscIndex3'; 'oscIndex4'; ...
                    'oscPeriod1Ms'; 'oscPeriod2Ms'; ...
                    'oscDurationSec'; ...
                    'nSpikesTotal'; 'nSpikesIn10s'; 'nSpikesInOsc'; ...
                    'nBurstsTotal'; 'nBurstsIn10s'; 'nBurstsInOsc'; ...
                    'nSpikesPerBurst'; 'nSpikesPerBurstIn10s'; ...
                    'nSpikesPerBurstInOsc'};
varLabelsAll = {'Oscillatory Index 1'; 'Oscillatory Index 2'; ...
                'Oscillatory Index 3'; 'Oscillatory Index 4'; ...
                'Oscillation Period 1 (ms)'; 'Oscillation Period 2 (ms)'; ...
                'Oscillation Duration (s)'; ...
                'Total Spike Count'; 'Number of Spikes in First 10 s'; 
                'Number of Spikes in Oscillation'; ...
                'Total Number of Bursts'; 'Number of Bursts in First 10 s'; ...
                'Number of Bursts in Oscillation'; ...
                'Number of Spikes Per Burst'; ...
                'Number of Spikes Per Burst in First 10 s'; ...
                'Number of Spikes Per Burst in Oscillation'};

%% Default values for optional arguments
plotTypeDefault = 'tuning';
computeChevronFlagDefault = false;
computeTimeTablesFlagDefault = false;
computePopAverageFlagDefault = false;
computeNormTablesFlagDefault = false;
computeNormChevronFlagDefault = false;
computeNormTimeTablesFlagDefault = false;
computeNormPopAverageFlagDefault = false;
computeSmoothTablesFlagDefault = false;
computeSmoothTimeTablesFlagDefault = false;
computeSmoothNormTimeTablesFlagDefault = false;
computeSmoothNormPopAverageFlagDefault = false;
plotAllFlagDefault = false;
plotByFileFlagDefault = [];
plotByPhaseFlagDefault = [];
plotPopAverageFlagDefault = [];
plotChevronFlagDefault = [];
plotNormByFileFlagDefault = [];
plotNormByPhaseFlagDefault = [];
plotNormPopAverageFlagDefault = [];
plotNormChevronFlagDefault = [];
plotSmoothByFileFlagDefault = [];
plotSmoothNormByFileFlagDefault = [];
plotSmoothNormPopAvgFlagDefault = [];
removeOutliersInPlotDefault = true;
directoryDefault = pwd;
inFolderDefault = '';                   % set later
outFolderDefault = '';                  % set later
figFolderDefault = '';                  % set later
phaseNumbersDefault = [];
sweepNumbersDefault = [];
sweepsRelToPhase2Default = [];
nSweepsLastOfPhaseDefault = 10;         % select from last 10 values by default
nSweepsToAverageDefault = 5;            % select 5 values by default
selectionMethodDefault = 'maxRange2Mean';   
                                        % select using maxRange2Mean by default
maxRange2MeanDefault = 40;              % range is not more than 40% of mean 
                                        %   by default
sweepLengthSecDefault = 60;             % sweep length is 60 seconds by default
timeLabelDefault = 'Time';
phaseLabelDefault = 'Phase';
phaseStrsDefault = {'Baseline', 'Wash-on', 'Wash-out'};
varsToPlotDefault = varsToPlotAll;
varLabelsDefault = varLabelsAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'ComputeChevronFlag', computeChevronFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeTimeTablesFlag', computeTimeTablesFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputePopAverageFlag', computePopAverageFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeNormTablesFlag', computeNormTablesFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeNormChevronFlag', computeNormChevronFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeNormTimeTablesFlag', ...
    computeNormTimeTablesFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeNormPopAverageFlag', ...
    computeNormPopAverageFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeSmoothTablesFlag', ...
    computeSmoothTablesFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeSmoothTimeTablesFlag', ...
    computeSmoothTimeTablesFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeSmoothNormTimeTablesFlag', ...
    computeSmoothNormTimeTablesFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ComputeSmoothNormPopAverageFlag', ...
    computeSmoothNormPopAverageFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotAllFlag', plotAllFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotByFileFlag', plotByFileFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotByPhaseFlag', plotByPhaseFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotPopAverageFlag', plotPopAverageFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotChevronFlag', plotChevronFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotNormByFileFlag', plotNormByFileFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotNormByPhaseFlag', plotNormByPhaseFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotNormPopAverageFlag', plotNormPopAverageFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotNormChevronFlag', plotNormChevronFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSmoothByFileFlag', plotSmoothByFileFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSmoothNormByFileFlag', plotSmoothNormByFileFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSmoothNormPopAvgFlag', plotSmoothNormPopAvgFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveOutliersInPlot', removeOutliersInPlotDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigFolder', figFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PhaseNumbers', phaseNumbersDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'SweepNumbers', sweepNumbersDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'SweepsRelToPhase2', sweepsRelToPhase2Default, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSweepsLastOfPhase', nSweepsLastOfPhaseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'NSweepsToAverage', nSweepsToAverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'SweepLengthSec', sweepLengthSecDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'TimeLabel', timeLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PhaseLabel', phaseLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PhaseStrings', phaseStrsDefault, ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['PhaseStrings must be a cell array of character arrays ', ...
                'or a string array!']));
addParameter(iP, 'VarsToPlot', varsToPlotDefault, ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['VarsToPlot must be a cell array of character arrays ', ...
                'or a string array!']));
addParameter(iP, 'VarLabels', varLabelsDefault, ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['VarLabels must be a cell array of character arrays ', ...
                'or a string array!']));

% Read from the Input Parser
parse(iP, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
computeChevronFlag = iP.Results.ComputeChevronFlag;
computeTimeTablesFlag = iP.Results.ComputeTimeTablesFlag;
computePopAverageFlag = iP.Results.ComputePopAverageFlag;
computeNormTablesFlag = iP.Results.ComputeNormTablesFlag;
computeNormChevronFlag = iP.Results.ComputeNormChevronFlag;
computeNormTimeTablesFlag = iP.Results.ComputeNormTimeTablesFlag;
computeNormPopAverageFlag = iP.Results.ComputeNormPopAverageFlag;
computeSmoothTablesFlag = iP.Results.ComputeSmoothTablesFlag;
computeSmoothTimeTablesFlag = iP.Results.ComputeSmoothTimeTablesFlag;
computeSmoothNormTimeTablesFlag = iP.Results.ComputeSmoothNormTimeTablesFlag;
computeSmoothNormPopAverageFlag = iP.Results.ComputeSmoothNormPopAverageFlag;
plotAllFlag = iP.Results.PlotAllFlag;
plotByFileFlag = iP.Results.PlotByFileFlag;
plotByPhaseFlag = iP.Results.PlotByPhaseFlag;
plotPopAverageFlag = iP.Results.PlotPopAverageFlag;
plotChevronFlag = iP.Results.PlotChevronFlag;
plotNormByFileFlag = iP.Results.PlotNormByFileFlag;
plotNormByPhaseFlag = iP.Results.PlotNormByPhaseFlag;
plotNormPopAverageFlag = iP.Results.PlotNormPopAverageFlag;
plotNormChevronFlag = iP.Results.PlotNormChevronFlag;
plotSmoothByFileFlag = iP.Results.PlotSmoothByFileFlag;
plotSmoothNormByFileFlag = iP.Results.PlotSmoothNormByFileFlag;
plotSmoothNormPopAvgFlag = iP.Results.PlotSmoothNormPopAvgFlag;
removeOutliersInPlot = iP.Results.RemoveOutliersInPlot;
directory = iP.Results.Directory;
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
figFolder = iP.Results.FigFolder;
phaseNumbers = iP.Results.PhaseNumbers;
sweepNumbers = iP.Results.SweepNumbers;
sweepsRelToPhase2 = iP.Results.SweepsRelToPhase2;
nSweepsLastOfPhase = iP.Results.NSweepsLastOfPhase;
nSweepsToAverage = iP.Results.NSweepsToAverage;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
maxRange2Mean = iP.Results.MaxRange2Mean;
sweepLengthSec = iP.Results.SweepLengthSec;
timeLabel = iP.Results.TimeLabel;
phaseLabel = iP.Results.PhaseLabel;
phaseStrs = iP.Results.PhaseStrings;
varsToPlot = iP.Results.VarsToPlot;
varLabels = iP.Results.VarLabels;

% Keep unmatched arguments for the plot_table() and plot_tuning_curve() function
otherArguments = iP.Unmatched;

% Force as a column cell array
[varsToPlot, varLabels] = argfun(@force_column_cell, varsToPlot, varLabels);

% Make sure the arguments match
if numel(varsToPlot) ~= numel(varLabels)
    error('VarsToPlot and VarLabels must have the same number of elements!');
end

%% Preparation
% Set default flags
fprintf('Setting default flags ...\n');
[plotByFileFlag, plotByPhaseFlag, plotPopAverageFlag, plotChevronFlag, ...
plotNormByFileFlag, plotNormByPhaseFlag, ...
plotNormPopAverageFlag, plotNormChevronFlag, ...
plotSmoothByFileFlag, plotSmoothNormByFileFlag, ...
plotSmoothNormPopAvgFlag] = ...
    argfun(@(x) set_default_flag(x, plotAllFlag), ...
                plotByFileFlag, plotByPhaseFlag, ...
                plotPopAverageFlag, plotChevronFlag, ...
                plotNormByFileFlag, plotNormByPhaseFlag, ...
                plotNormPopAverageFlag, plotNormChevronFlag, ...
                plotSmoothByFileFlag, plotSmoothNormByFileFlag, ...
                plotSmoothNormPopAvgFlag);

% Set compute flags
if plotPopAverageFlag
    computePopAverageFlag = true;
end
if plotNormPopAverageFlag
    computeNormPopAverageFlag = true;
end
if plotSmoothNormPopAvgFlag
    computeSmoothNormPopAverageFlag = true;
end
if plotByFileFlag || plotByPhaseFlag || computePopAverageFlag
    computeTimeTablesFlag = true;
end
if plotNormByFileFlag || plotNormByPhaseFlag || computeNormPopAverageFlag
    computeNormTimeTablesFlag = true;
end
if plotSmoothNormByFileFlag || computeSmoothNormPopAverageFlag
    computeSmoothNormTimeTablesFlag = true;
end
if plotNormChevronFlag
    computeNormChevronFlag = true;
end
if computeNormTimeTablesFlag
    % TODO: redundant
    computeNormTablesFlag = true;
end
if plotSmoothByFileFlag || computeSmoothNormTimeTablesFlag
    computeSmoothTimeTablesFlag = true;
end
if plotSmoothByFileFlag || computeSmoothTimeTablesFlag
    computeSmoothTablesFlag = true;
end
if plotChevronFlag || computeNormChevronFlag || ...
        computeNormTablesFlag || computeSmoothNormTimeTablesFlag
    computeChevronFlag = true;
end

% Decide on the input directory
if isempty(inFolder)
    inFolder = directory;
end

% Decide on the output directory
if isempty(outFolder)
    outFolder = fullfile(inFolder, [create_time_stamp, '_', outDirSuffix]);
end

% Decide on the figure directory
if isempty(figFolder)
    figFolder = outFolder;
end

% Find all files with the pattern *slice*_params in the file name
fprintf('Finding all spreadsheets ...\n');
[~, sliceParamSheets] = all_files('Directory', inFolder, 'Keyword', 'slice', ...
                                'Suffix', 'params', 'Extension', 'csv', ...
                                'ForceCellOutput', true);

% If nothing found, return
if isempty(sliceParamSheets)
    fprintf('There are no slice-dependent spreadsheets in %s!\n', inFolder);
    return
end

% Extract the common prefix
prefix = extract_fileparts(sliceParamSheets, 'commonprefix');

% If no common prefix, use the directory name as the prefix
if isempty(prefix)
    prefix = extract_common_directory(sliceParamSheets, 'BaseNameOnly', true);
end

% Modify prefix if necessary
if ~isempty(phaseNumbers)
    % Concatenate phase numbers into a string
    phaseNumbersString = num2str(phaseNumbers, '%d');

    % Append the phase numbers to the prefix
    prefix = [prefix, '_phases', phaseNumbersString];
end
if ~isempty(sweepNumbers)
    % Create a sweep number string
    sweepNumbersString = create_label_from_sequence(sweepNumbers);

    % Append the phase numbers to the prefix
    prefix = [prefix, '_sweeps', sweepNumbersString];
end
if ~isempty(sweepsRelToPhase2)
    % Create a sweep number string
    sweepsRelToPhase2String = create_label_from_sequence(sweepsRelToPhase2);

    % Append the phase numbers to the prefix
    prefix = [prefix, '_sweeps', sweepsRelToPhase2String, 'relToPhase2'];
end

% Extract the distinct parts of the file names
fileLabels = extract_fileparts(sliceParamSheets, 'distinct');

% Generate variable labels for normalized plots
if plotNormByFileFlag || plotNormByPhaseFlag || ...
        plotNormChevronFlag || plotNormPopAverageFlag || ...
        plotSmoothNormByFileFlag || plotSmoothNormPopAvgFlag
    varLabelsNorm = repmat({'% of baseline'}, size(varLabels));
end

% Create table labels
tableLabels = strcat(prefix, {': '}, varLabels);

% Create table names
tableNames = strcat(prefix, '_', varsToPlot);

% Create figure names
[figNamesAvgd, figNamesByFile, figNamesByPhase, figNamesPopAvg, ...
    figNamesNormAvgd, figNamesNormByFile, ...
    figNamesNormByPhase, figNamesNormPopAvg, ...
    figNamesSmoothByFile, figNamesSmoothNormByFile, ...
    figNamesSmoothNormPopAvg] = ...
    argfun(@(x) fullfile(figFolder, strcat(tableNames, '_', x, '.epsc')), ...
            'chevron', 'byFile', 'byPhase', 'popAverage', ...
            'normChevron', 'normByFile', 'normByPhase', 'normPopAverage', ...
            'smoothByFile', 'smoothNormByFile', ...
            'smoothNormPopAverage');

% Create paths for spreadsheet files
[chevronTablePaths, popAvgTablePaths, smoothTablePaths, ...
    normTablePaths, normChevronTablePaths, normPopAvgTablePaths, ...
    smoothNormTablePaths, smoothNormPopAvgTablePaths] = ...
    argfun(@(x) fullfile(outFolder, strcat(tableNames, '_', x, '.csv')), ...
            'chevron', 'popAverage', 'smooth', ...
            'normalized', 'normChevron', 'normPopAverage', ...
            'smoothNorm', 'smoothNormPopAverage');

% Create paths for mat files
[measureTablesMatPath, measureTimeTablesMatPath, chevronTablesMatPath, ...
    popAvgTablesMatPath, normTimeTablesMatPath, ...
    normalizedChevronTablesMatPath, normalizedPopAvgTablesMatPath, ...
    smoothTablesMatPath, smoothTimeTablesMatPath, ...
    smoothNormTimeTablesMatPath, smoothNormPopAvgTablesMatPath] = ...
    argfun(@(x) fullfile(outFolder, [prefix, '_', x, '.mat']), ...
            'measureTables', 'measureTimeTables', 'chevronTables', ...
            'popAvgTables', 'normalizedTimeTables', ...
            'normalizedChevronTables', 'normalizedPopAvgTables', ...
            'smoothTables', 'smoothTimeTables', ...
            'smoothNormTimeTables', 'smoothNormPopAvgTables');

% Check if output directories exist
check_dir({outFolder, figFolder});

%% Read in data and preprocess
% Read all slice parameter spreadsheets
fprintf('Reading measure spreadsheets ...\n');
sliceParamsTables = cellfun(@readtable, sliceParamSheets, ...
                            'UniformOutput', false);

% Backup the spreadsheets in the figure folder
if backupSheets
    copy_into(sliceParamSheets, figFolder);
end

% Create a time column (time in minutes since drug on)
fprintf('Creating time columns ...\n');
sliceParamsTables = ...
    cellfun(@(x) create_time_rel_to_drugon(x, sweepLengthSec), ...
                sliceParamsTables, 'UniformOutput', false);

% Combine with phase number information
%   Note: each row results in a different table
varsToCombine = [varsToPlot, repmat({'phaseNumber'}, size(varsToPlot))];

% Create the phaseNumber variables for the combined tables
phaseVars = strcat('phaseNumber_', fileLabels);

% Restrict to certain phases if requested
if ~isempty(phaseNumbers)
    fprintf('Restricting to phases %s ...\n', num2str(phaseNumbers, '%d, '));
    sliceParamsTables = ...
        cellfun(@(x) x(ismember(x.phaseNumber, phaseNumbers), :), ...
                sliceParamsTables, 'UniformOutput', false);
end
if ~isempty(sweepNumbers)
    fprintf('Restricting to sweeps %s ...\n', sweepNumbersString);
    sliceParamsTables = ...
        cellfun(@(x) x(ismember(x.sweepNumber, sweepNumbers), :), ...
                sliceParamsTables, 'UniformOutput', false);
end
if ~isempty(sweepsRelToPhase2)
    fprintf('Restricting to sweeps %s relative to phase 2 ...\n', ...
            sweepsRelToPhase2String);

    sliceParamsTables = ... 
        cellfun(@(x) restrict_to_relative_sweeps(x, sweepsRelToPhase2, 2), ...
                sliceParamsTables, 'UniformOutput', false);
end

% Generate phase tick locations
if isempty(phaseNumbers)
    phaseTickLocs = 1:3;
else
    phaseTickLocs = 1:numel(phaseNumbers);
    phaseStrs = phaseStrs(phaseNumbers);
end

%% Combine variables across tables
fprintf('Combining variables across tables ...\n');
measureTables = combine_variables_across_tables(sliceParamsTables, ...
                'Keys', 'Time', 'VariableNames', varsToCombine, ...
                'InputNames', fileLabels, 'OmitVarName', false, ...
                'OutFolder', outFolder, 'Prefix', prefix, 'SaveFlag', true);

save(measureTablesMatPath, 'measureTables', '-mat');

%% Compute smoothed tables
if computeSmoothTablesFlag
    smoothTables = ...
        cellfun(@(x, y, z) modify_table(x, smoothFunc, ...
                                        'VariableNames', y, 'SheetName', z), ...
                measureTables, varsToPlot, smoothTablePaths, ...
                'UniformOutput', false);

    save(smoothTablesMatPath, 'smoothTables', '-mat');
end

%% Average over each phase
if computeChevronFlag
    % Average over the last nSweepsToAverage sweeps
    fprintf('Creating Chevron tables ...\n');
    chevronTables = ...
        cellfun(@(x, y) average_last_of_each_phase(x, nSweepsLastOfPhase, ...
                                    nSweepsToAverage, selectionMethod, ...
                                    maxRange2Mean, y), ...
                        measureTables, chevronTablePaths, 'UniformOutput', false);

    save(chevronTablesMatPath, 'chevronTables', '-mat');
end

%% Normalize raw data to baseline
% TODO: Don't need this. Normalize time tables directly?
if computeNormTablesFlag
    fprintf('Computing normalized tables ...\n');
    normalizedMeasureTables = ...
        cellfun(@(x, y, z, w) normalize_to_baseline(x, y, z, w), ...
                measureTables, chevronTables, varsToPlot, normTablePaths, ...
                'UniformOutput', false);
end

%% Normalize averaged data to baseline
if computeNormChevronFlag
    fprintf('Normalizing to baseline ...\n');
    normalizedChevronTables = ...
        cellfun(@(x, y) normalize_to_first_row(x, y), ...
                    chevronTables, normChevronTablePaths, 'UniformOutput', false);

    save(normalizedChevronTablesMatPath, 'normalizedChevronTables', '-mat');
end

%% Convert to timetables
if computeTimeTablesFlag
    fprintf('Converting measure tables to time tables ...\n');
    measureTimeTables = cellfun(@table2timetable, ...
                                measureTables, 'UniformOutput', false);

    save(measureTimeTablesMatPath, 'measureTimeTables', '-mat');
end

%% Convert normalized data to timetables
if computeNormTimeTablesFlag
    fprintf('Converting normalized measure tables to time tables ...\n');
    normTimeTables = cellfun(@table2timetable, ...
                            normalizedMeasureTables, 'UniformOutput', false);
    save(normTimeTablesMatPath, 'normTimeTables', '-mat');
end

%% Convert smoothed data to timetables
if computeSmoothTimeTablesFlag
    fprintf('Converting smoothed measure tables to time tables ...\n');
    smoothTimeTables = cellfun(@table2timetable, ...
                            smoothTables, 'UniformOutput', false);
    save(smoothTimeTablesMatPath, 'smoothTimeTables', '-mat');
end

%% Normalize smoothed data to baseline
if computeSmoothNormTimeTablesFlag
    fprintf('Computing smoothed, then normalized tables ...\n');
    smoothNormTimeTables = ...
        cellfun(@(x, y, z, w) normalize_to_baseline(x, y, z, w), ...
                smoothTimeTables, chevronTables, ...
                varsToPlot, smoothNormTablePaths, ...
                'UniformOutput', false);
end

%% Average over slices
if computePopAverageFlag
    fprintf('Computing population averages ...\n');
    popAvgTables = cellfun(@(x, y, z) compute_population_average(x, ...
                                            'VarStr', y, 'SheetName', z), ...
                            measureTimeTables, varsToPlot, popAvgTablePaths, ...
                            'UniformOutput', false);

    save(popAvgTablesMatPath, 'popAvgTables', '-mat');
end

%% Average over slices for normalized data
if computeNormPopAverageFlag
    fprintf('Computing population averages of normalized data ...\n');
    normPopAvgTables = cellfun(@(x, y, z) compute_population_average(x, ...
                                            'VarStr', y, 'SheetName', z), ...
                            normTimeTables, varsToPlot, normPopAvgTablePaths, ...
                            'UniformOutput', false);

    save(normalizedPopAvgTablesMatPath, 'normPopAvgTables', '-mat');
end

%% Average over slices for smoothed, then normalized data
if computeSmoothNormPopAverageFlag
    fprintf('Computing population averages of smoothed, then normalized data ...\n');
    smoothNormPopAvgTables = cellfun(@(x, y, z) compute_population_average(x, ...
                                            'VarStr', y, 'SheetName', z), ...
                            smoothNormTimeTables, varsToPlot, ...
                            smoothNormPopAvgTablePaths, ...
                            'UniformOutput', false);

    save(smoothNormPopAvgTablesMatPath, 'smoothNormPopAvgTables', '-mat');
end

%% Plot Chevron plots
if plotChevronFlag || plotNormChevronFlag
    % Generate figure titles for Chevron plots
    figTitlesChevron = strcat(varLabels, [' avg over ', ...
                                num2str(nSweepsToAverage), ' of last ', ...
                                num2str(nSweepsLastOfPhase), ' sweeps']);
end

% Plot regular Chevron plots
if plotChevronFlag
    % Generate variable labels
    varLabelsChevron = varLabels;

    close all;
    fprintf('Plotting Chevron plots ...\n');
    cellfun(@(a, b, c, d, e, f) plot_chevron_from_table(a, b, c, d, e, f, ...
                                    plotType, fileLabels, ...
                                    phaseTickLocs, phaseStrs, phaseLabel, ...
                                    plotForCorel, otherArguments), ...
            chevronTables, varsToPlot, varLabelsChevron, ...
            tableLabels, figTitlesChevron, figNamesAvgd);
end

% Plot Normalized Chevron plots
if plotNormChevronFlag
    close all;
    fprintf('Plotting Normalized Chevron plots ...\n');
    cellfun(@(x, y, z, w, v, u) plot_table(x, 'PlotMode', 'overlapped', ...
                                'PlotType', plotType, ...
                                'VarsToPlot', strcat(y, '_', fileLabels), ...
                                'PTicks', phaseTickLocs, ...
                                'PTickLabels', phaseStrs, ...
                                'ReadoutLabel', z, 'TableLabel', w, ...
                                'PLabel', phaseLabel, ...
                                'FigTitle', v, 'FigName', u, ...
                                'LegendLocation', 'eastoutside', ...
                                'RemoveOutliers', false, ...
                                'Marker', 'o', 'MarkerFaceColor', 'auto', ...
                                otherArguments), ...
            normalizedChevronTables, varsToPlot, varLabelsNorm, ...
            tableLabels, figTitlesChevron, figNamesNormAvgd);
end

%% Plot each column with a different color
if plotByFileFlag
    close all;
    fprintf('Plotting each column with a different color ...\n');
    cellfun(@(x, y, z, w, v) plot_table(x, 'PlotMode', 'overlapped', ...
                'PhaseVariables', phaseVars, 'PhaseLabels', phaseStrs, ...
                'PlotPhaseBoundaries', true, 'FigExpansion', 2, ...
                'PlotPhaseAverages', true, 'PlotIndSelected', true, ...
                'PlotAverageWindows', false, ...
                'NLastOfPhase', nSweepsLastOfPhase, ...
                'NToAverage', nSweepsToAverage, ...
                'SelectionMethod', selectionMethod, ...
                'MaxRange2Mean', maxRange2Mean, ...
                'ColorByPhase', false, ...
                'PlotType', plotType, ...
                'VarsToPlot', strcat(y, '_', fileLabels), ...
                'ReadoutLabel', z, 'TableLabel', w, ...
                'PLabel', timeLabel, 'FigName', v, ...
                'RemoveOutliers', removeOutliersInPlot, otherArguments), ...
            measureTimeTables, varsToPlot, varLabels, tableLabels, figNamesByFile);
end

if plotNormByFileFlag
    close all;
    fprintf('Plotting each column with a different color, data normalized ...\n');
    cellfun(@(x, y, z, w, v) plot_table(x, 'PlotMode', 'overlapped', ...
                'PhaseVariables', phaseVars, 'PhaseLabels', phaseStrs, ...
                'PlotPhaseBoundaries', true, 'PlotPhaseAverages', false, ...
                'PlotIndSelected', false, 'ColorByPhase', false, ...
                'PlotAverageWindows', false, ...
                'NLastOfPhase', nSweepsLastOfPhase, ...
                'NToAverage', nSweepsToAverage, ...
                'SelectionMethod', selectionMethod, ...
                'MaxRange2Mean', maxRange2Mean, ...
                'PlotType', plotType, ...
                'VarsToPlot', strcat(y, '_', fileLabels), ...
                'ReadoutLabel', z, 'TableLabel', w, ...
                'PLabel', timeLabel, 'FigName', v, ...
                'RemoveOutliers', removeOutliersInPlot, otherArguments), ...
            normTimeTables, varsToPlot, varLabelsNorm, tableLabels, figNamesNormByFile);
end

if plotSmoothByFileFlag
    close all;
    tableLabelsSmooth = strcat(tableLabels, ' smoothed');

    fprintf('Plotting each column with a different color, data smoothed ...\n');
    cellfun(@(x, y, z, w, v) plot_table(x, 'PlotMode', 'overlapped', ...
                'PhaseVariables', phaseVars, 'PhaseLabels', phaseStrs, ...
                'PlotPhaseBoundaries', true, 'PlotPhaseAverages', false, ...
                'PlotIndSelected', false, 'ColorByPhase', false, ...
                'PlotAverageWindows', false, ...
                'NLastOfPhase', nSweepsLastOfPhase, ...
                'NToAverage', nSweepsToAverage, ...
                'SelectionMethod', selectionMethod, ...
                'MaxRange2Mean', maxRange2Mean, ...
                'PlotType', plotType, ...
                'VarsToPlot', strcat(y, '_', fileLabels), ...
                'ReadoutLabel', z, 'TableLabel', w, ...
                'PLabel', timeLabel, 'FigName', v, ...
                'RemoveOutliers', removeOutliersInPlot, otherArguments), ...
            smoothTimeTables, varsToPlot, varLabels, ...
            tableLabelsSmooth, figNamesSmoothByFile);
end

if plotSmoothNormByFileFlag
    close all;
    tableLabelsSmoothNorm = strcat(tableLabels, ' smoothed, then normalized');

    fprintf(['Plotting each column with a different color, ', ...
                'data smoothed then normalized ...\n']);
    cellfun(@(x, y, z, w, v) plot_table(x, 'PlotMode', 'overlapped', ...
                'PhaseVariables', phaseVars, 'PhaseLabels', phaseStrs, ...
                'PlotPhaseBoundaries', true, 'PlotPhaseAverages', false, ...
                'PlotIndSelected', false, 'ColorByPhase', false, ...
                'PlotAverageWindows', false, ...
                'NLastOfPhase', nSweepsLastOfPhase, ...
                'NToAverage', nSweepsToAverage, ...
                'SelectionMethod', selectionMethod, ...
                'MaxRange2Mean', maxRange2Mean, ...
                'PlotType', plotType, ...
                'VarsToPlot', strcat(y, '_', fileLabels), ...
                'ReadoutLabel', z, 'TableLabel', w, ...
                'PLabel', timeLabel, 'FigName', v, ...
                'RemoveOutliers', removeOutliersInPlot, otherArguments), ...
            smoothNormTimeTables, varsToPlot, varLabelsNorm, ...
            tableLabelsSmoothNorm, figNamesSmoothNormByFile);
end

%% Plot each phase with a different color
if plotByPhaseFlag
    close all;
    fprintf('Plotting each phase with a different color ...\n');
    cellfun(@(x, y, z, w, v) plot_table(x, 'PlotMode', 'overlapped', ...
                'PhaseVariables', phaseVars, 'PhaseLabels', phaseStrs, ...
                'PlotPhaseBoundaries', false, 'FigExpansion', 2, ...
                'PlotPhaseAverages', true, 'PlotIndSelected', true, ...
                'PlotAverageWindows', false, ...
                'NLastOfPhase', nSweepsLastOfPhase, ...
                'NToAverage', nSweepsToAverage, ...
                'SelectionMethod', selectionMethod, ...
                'MaxRange2Mean', maxRange2Mean, ...
                'ColorByPhase', true, ...
                'PlotType', plotType, ...
                'VarsToPlot', strcat(y, '_', fileLabels), ...
                'ReadoutLimits', [0, Inf], ...
                'ReadoutLabel', z, 'TableLabel', w, ...
                'PLabel', timeLabel, 'FigName', v, ...
                'RemoveOutliers', removeOutliersInPlot, otherArguments), ...
            measureTimeTables, varsToPlot, varLabels, tableLabels, figNamesByPhase);

    %                                 'PhaseLabels', phaseStrs, ...
end

if plotNormByPhaseFlag
    close all;
    fprintf('Plotting each phase with a different color, data normalized ...\n');
    cellfun(@(x, y, z, w, v) plot_table(x, 'PlotMode', 'overlapped', ...
                'PhaseVariables', phaseVars, 'PhaseLabels', phaseStrs, ...
                'PlotPhaseBoundaries', false, 'PlotPhaseAverages', false, ...
                'PlotIndSelected', false, 'ColorByPhase', true, ...
                'PlotAverageWindows', false, ...
                'NLastOfPhase', nSweepsLastOfPhase, ...
                'NToAverage', nSweepsToAverage, ...
                'SelectionMethod', selectionMethod, ...
                'MaxRange2Mean', maxRange2Mean, ...
                'PlotType', plotType, ...
                'VarsToPlot', strcat(y, '_', fileLabels), ...
                'ReadoutLimits', [0, Inf], ...
                'ReadoutLabel', z, 'TableLabel', w, ...
                'PLabel', timeLabel, 'FigName', v, ...
                'RemoveOutliers', removeOutliersInPlot, otherArguments), ...
            normTimeTables, varsToPlot, varLabelsNorm, tableLabels, figNamesNormByPhase);
end

%% Plot population averages
if plotPopAverageFlag
    % Create titles
    figTitlesPopAvg = replace(tableLabels, '_', '\_');

    close all;
    fprintf('Plotting population averages ...\n');
    cellfun(@(x, y, z, w, v) ...
                plot_tuning_curve(x.Properties.RowTimes, x{:, [y, '_mean']}, ...
                    'LowerCI', x{:, [y, '_lower95']}, ...
                    'UpperCI', x{:, [y, '_upper95']}, ...
                    'PhaseVectors', x{:, 'phaseNumber'}, ...
                    'PhaseLabels', phaseStrs, ...
                    'PBoundaryType', 'verticalShade', 'ColorByPhase', false, ...
                    'PlotPhaseBoundaries', true, 'PlotPhaseAverages', false, ...
                    'PlotIndSelected', false, 'PlotAverageWindows', true, ...
                    'ClearFigure', true, ...
                    'ReadoutLimits', [0, Inf], ...
                    'PLabel', timeLabel, 'ReadoutLabel', z, ...
                    'FigTitle', w, 'FigName', v, otherArguments), ...
            popAvgTables, varsToPlot, varLabels, ...
            figTitlesPopAvg, figNamesPopAvg);

end

%% Plot population averages of normalized data
if plotNormPopAverageFlag
    % Create titles
    figTitlesPopAvg = replace(tableLabels, '_', '\_');

    close all;
    fprintf('Plotting population averages of normalized data ...\n');
    cellfun(@(x, y, z, w, v) ...
                plot_tuning_curve(x.Properties.RowTimes, x{:, [y, '_mean']}, ...
                    'LowerCI', x{:, [y, '_lower95']}, ...
                    'UpperCI', x{:, [y, '_upper95']}, ...
                    'PhaseVectors', x{:, 'phaseNumber'}, ...
                    'PhaseLabels', phaseStrs, ...
                    'PBoundaryType', 'verticalShade', 'ColorByPhase', false, ...
                    'PlotPhaseBoundaries', true, 'PlotPhaseAverages', false, ...
                    'PlotIndSelected', false, 'PlotAverageWindows', true, ...
                    'ClearFigure', true, ...
                    'ReadoutLimits', [0, Inf], ...
                    'PLabel', timeLabel, 'ReadoutLabel', z, ...
                    'FigTitle', w, 'FigName', v, otherArguments), ...
            normPopAvgTables, varsToPlot, varLabelsNorm, ...
            figTitlesPopAvg, figNamesNormPopAvg);
end

%% Plot population averages of smoothed, then normalized data 
if plotSmoothNormPopAvgFlag
    % Create titles
    figTitlesPopAvg = replace(tableLabels, '_', '\_');

    close all;
    fprintf('Plotting population averages of smoothed, then normalized data  ...\n');
    cellfun(@(x, y, z, w, v) ...
                plot_tuning_curve(x.Properties.RowTimes, x{:, [y, '_mean']}, ...
                    'LowerCI', x{:, [y, '_lower95']}, ...
                    'UpperCI', x{:, [y, '_upper95']}, ...
                    'PhaseVectors', x{:, 'phaseNumber'}, ...
                    'PhaseLabels', phaseStrs, ...
                    'PBoundaryType', 'verticalShade', 'ColorByPhase', false, ...
                    'PlotPhaseBoundaries', true, 'PlotPhaseAverages', false, ...
                    'PlotIndSelected', false, 'PlotAverageWindows', true, ...
                    'ClearFigure', true, ...
                    'ReadoutLimits', [0, Inf], ...
                    'PLabel', timeLabel, 'ReadoutLabel', z, ...
                    'FigTitle', w, 'FigName', v, 'FigTypes', figTypes, ...
                    otherArguments), ...
            smoothNormPopAvgTables, varsToPlot, varLabelsNorm, ...
            figTitlesPopAvg, figNamesSmoothNormPopAvg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_chevron_from_table(chevronTable, varToPlot, varLabel, ...
                                    tableLabel, figTitle, figName, ...
                                    plotType, fileLabels, ...
                                    phaseTickLocs, phaseStrs, phaseLabel, ...
                                    plotForCorel, otherArguments)

%% Hard-coded parameters
chevronFigWidth = 3.25;           % in cm
chevronFigHeight = 4;             % in cm
chevronMarkerSize = 1;          % marker size in points
barInsetWidth = 1;              % figure width in cm
barInsetHeight = 2;             % figure height in cm
barInsetRulerLineWidth = 0.5;
barInsetAxisFontSize = 5;
barInsetLabelsFontSize = 6;
figTypes = {'png', 'epsc'};

%% Preparation
figPathBase = extract_fileparts(figName, 'pathbase');

%% Plot
% Create a new figure
fig = set_figure_properties('AlwaysNew', true);

% Plot the Chevron plot
%{
plot_table(chevronTable, 'PlotMode', 'overlapped', 'PlotType', plotType, ...
            'VarsToPlot', strcat(varToPlot, '_', fileLabels), ...
            'PTicks', phaseTickLocs, 'PTickLabels', phaseStrs, ...
            'ReadoutLabel', varLabel, 'TableLabel', tableLabel, ...
            'PLabel', phaseLabel, 'FigTitle', figTitle, ...
            'LegendLocation', 'eastoutside', 'RemoveOutliers', false, ...
            'RunTTest', true, 'RunRankTest', true, ...
            'Marker', 'o', 'MarkerFaceColor', 'auto', ...
            otherArguments);
%}

varsToPlot = strcat(varToPlot, '_', fileLabels);
chevronData = chevronTable(:, varsToPlot);

plot_chevron(chevronData, 'PlotMeanValues', true, ...
            'PlotMeanDifference', true, 'PlotErrorBars', false, ...
            'PTicks', phaseTickLocs, 'PTickLabels', phaseStrs, ...
            'ReadoutLabel', varLabel, 'PLabel', phaseLabel, ...
            'FigTitle', figTitle, 'LegendLocation', 'eastoutside', ...
            'RunTTest', true, 'RunRankTest', true, ...
            'Marker', 'o', 'MarkerFaceColor', 'auto', ...
            otherArguments);

% Update for CorelDraw
% TODO: Move this into plot_chevron.m?
if plotForCorel
    % Create path for original figure
    figNameOrig = [figPathBase, '_orig'];

    % Save original figure as png
    save_all_figtypes(fig, figNameOrig, 'png');

    % Update figure for CorelDraw
    update_figure_for_corel(fig, 'Units', 'centimeters', ...
                    'Width', chevronFigWidth, 'Height', chevronFigHeight, ...
                    'PlotMarkerSize', chevronMarkerSize, ...
                    'RemoveLegends', true);
end

% Save figure
save_all_figtypes(fig, figName, figTypes);

% Create figure for bar inset
fig2 = set_figure_properties('AlwaysNew', true);

% Create figure path base for bar inset
figPathBaseBar = [figPathBase, '_bar_inset'];

% Plot bar inset
plot_chevron_bar_inset(chevronTable, 'PTickLabels', phaseStrs);

% Update for CorelDraw
% TODO: Move this into plot_chevron_bar_inset.m?
if plotForCorel
    % Create path for original figure
    figPathBaseBarOrig = [figPathBase, '_bar_inset_orig'];

    % Save figure
    save_all_figtypes(fig2, figPathBaseBarOrig, 'png');

    % Update figure for CorelDraw
    update_figure_for_corel(fig2, 'Units', 'centimeters', ...
                        'Width', barInsetWidth, 'Height', barInsetHeight, ...
                        'RulerLineWidth', barInsetRulerLineWidth, ...
                        'AxisFontSize', barInsetAxisFontSize, ...
                        'LabelsFontSize', barInsetLabelsFontSize);
end

% Save figure
save_all_figtypes(fig2, figPathBaseBar, figTypes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myTable = create_time_rel_to_drugon(myTable, sweepLengthSec)
%% Creates a time column in minutes since drug onset
% Count the number of rows
nRows = height(myTable);

% Get the phase numbers
phaseNum = myTable.phaseNumber;

% Get the first row that is drug on
%   Note: this is set #2
rowDrugOn = find(phaseNum == 2, 1, 'first');

% Error check
if isempty(rowDrugOn)
    error('Phase 2 not found in one of the slices!');
end

% Compute time in sweeps relative to drug on
timeSwps = transpose(1:nRows) - rowDrugOn;

% Convert to minutes
timeMin = timeSwps / (sweepLengthSec / 60);

% Convert to a duration vector
Time = minutes(timeMin);

% Add a time column to the table
myTable = addvars(myTable, Time, 'Before', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outTable = ...
            average_last_of_each_phase(inTable, nSweepsLastOfPhase, ...
                                nSweepsToAverage, selectionMethod, ...
                                maxRange2Mean, sheetPath)
%% Averages over the last nSweepsToAverage sweeps of each phase
% Note: all distinct identifiers must have a matching phase variable column
% TODO: Fix

% Remove the Time column if exists
if is_var_in_table('Time', inTable)
    inTable = removevars(inTable, 'Time');
end

% Get all variable names that are left
allVarNames = transpose(inTable.Properties.VariableNames);

% Find the phase variable names
[~, phaseVars] = find_in_strings('phase.*', allVarNames, 'SearchMode', 'reg');

% Collect the rest of the variable names
readoutVars = setdiff(allVarNames, phaseVars);

% Extract the unique identifiers
uniqueIds = extract_distinct_fileparts(readoutVars, 'Delimiter', '_');

% Find matching phase variables for each unique identifiers
[~, phaseVars] = ...
    cellfun(@(x) find_in_strings(x, phaseVars, 'SearchMode', 'substrings'), ...
            uniqueIds, 'UniformOutput', false);

% Find the row indices for the last nSweepsLastOfPhase sweeps for each phase,
%   for each unique phase ID
indLastOfPhase = cellfun(@(x) find_last_ind_each_phase(inTable{:, x}, ...
                                                nSweepsLastOfPhase), ...
                    phaseVars, 'UniformOutput', false);

% Count the number of phases for each unique phase ID
nPhases = count_vectors(indLastOfPhase);

% Compute the maximum number of phases
maxNPhases = max(nPhases);

% Creat a phase number column
phaseNumber = transpose(1:maxNPhases);

% Average the nSweepsToAverage sweeps within maxRange2Mean in
%    the last nSweepsLastOfPhase sweeps of each phase
readoutAvg = ...
    cellfun(@(x, y) cellfun(@(z) compute_phase_average(inTable{:, x}, ...
                                        'ReturnLastTrial', true, ...
                                        'Indices', z, ...
                                        'NToAverage', nSweepsToAverage, ...
                                        'SelectionMethod', selectionMethod, ...
                                        'MaxRange2Mean', maxRange2Mean), ...
                            y), ...
            readoutVars, indLastOfPhase, 'UniformOutput', false);

% Force as a matrix, padding each vector with NaNs at the end if necessary
readoutAvg = force_matrix(readoutAvg, 'AlignMethod', 'leftAdjustPad');

% Create an averaged table
outTable = array2table(readoutAvg);
outTable = addvars(outTable, phaseNumber, 'Before', 1);
outTable.Properties.VariableNames = vertcat({'phaseNumber'}, readoutVars);

% Save the table
writetable(outTable, sheetPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indLastEachPhase = find_last_ind_each_phase(phaseVec, nLastIndices)
%% Find the last nLastIndices indices for each phase in the phaseVec

% Get the unique phases
uniquePhases = unique_custom(phaseVec, 'stable', 'IgnoreNaN', true);

% Find the last index for each phase
lastIndexEachPhase = arrayfun(@(x) find(ismatch(phaseVec, x), 1, 'last'), ...
                                uniquePhases);

% Compute the first index for each phase
firstIndexEachPhase = lastIndexEachPhase - nLastIndices + 1;

% Construct the last nLastIndices indices
indLastEachPhase = create_indices('IndexStart', firstIndexEachPhase, ...
                                    'IndexEnd', lastIndexEachPhase, ...
                                    'MaxNum', nLastIndices, ...
                                    'ForceCellOutput', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normalizedTable = normalize_to_baseline(table, avgTable, varStr, sheetPath)
%% Normalizes all values to baseline
% TODO: Pull out to its own function

% Get all the variable names
columnNames = table.Properties.VariableNames;

% Select the columns that contains given variable string
columnsToNormalize = find(contains(columnNames, varStr));

% Run through all columns containing given variable string
normalizedTable = table;
for idxCol = columnsToNormalize
    % Get the column name
    thisColumnName = columnNames{idxCol};

    % Extract the baseline value for this column
    baseValue = avgTable{1, thisColumnName};

    % Extract this column
    thisColumn = table{:, idxCol};
    
    % Normalize the column with the baseline value
    normalizedTable{:, idxCol} = (thisColumn ./ baseValue) .* 100;
end

% Save the table
if istimetable(normalizedTable)
    write_timetable(normalizedTable, sheetPath);
else
    writetable(normalizedTable, sheetPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normalizedTable = normalize_to_first_row(table, sheetPath)
%% Normalizes all values to the first row
% TODO: Pull out to its own function

% Extract the first row
firstRow = table{1, :};
% firstRow = table{3, :};

% Count the number of rows
nRows = height(table);

% Divide each row by the first row and multiply by 100 %
% normalizedTable = rowfun(@(x) x ./ firstRow, table);
for iRow = 1:nRows
    table{iRow, :} = (table{iRow, :} ./ firstRow) .* 100;
end

normalizedTable = table;

% Save the table
writetable(normalizedTable, sheetPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paramsTable = restrict_to_relative_sweeps (paramsTable, ...
                                        relativeSweepNums, phaseNumber)

% Extract from the parameters table
phaseNumbers = paramsTable.phaseNumber;
sweepNumbers = paramsTable.sweepNumber;

% Find the index of the last sweep from the previous phase
idxLastSweepPrevPhase = find(phaseNumbers < phaseNumber, 1, 'last');
if isempty(idxLastSweepPrevPhase)
    idxLastSweepPrevPhase = 0;
end

% Shift the relative sweep numbers to the original
sweepNumbersToRestrictTo = idxLastSweepPrevPhase + relativeSweepNums;

% Restrict to those sweeps
paramsTable = paramsTable(ismember(sweepNumbers, sweepNumbersToRestrictTo), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find the row indices for the last nSweepsToAverage sweeps for each phase,
%   for each unique phase ID
indToAvg = cellfun(@(x) find_last_ind_each_phase(inTable{:, x}, nSweepsToAverage), ...
                    phaseVars, 'UniformOutput', false);
% Average the last nSweepsToAverage sweeps for each phase
readoutAvg = cellfun(@(x, y) cellfun(@(z) nanmean(inTable{:, x}(z)), y), ...
                    readoutVars, indToAvg, 'UniformOutput', false);

% Count the number of unique phases
nPhases = numel(uniquePhases);

phaseStrs = {'baseline', 'washon', 'washoff'};

% File patterns
sliceFilePattern = '.*slice.*';

chevronTablePaths = fullfile(outFolder, strcat(tableNames, '_chevron.csv'));
normChevronTablePaths = fullfile(outFolder, strcat(tableNames, '_normChevron.csv'));
popAvgTablePaths = fullfile(outFolder, strcat(tableNames, '_popAverage.csv'));

measureTablesMatPath = fullfile(outFolder, [prefix, '_', 'measureTables.mat']);
chevronTablesMatPath = fullfile(outFolder, [prefix, '_', 'chevronTables.mat']);
normalizedChevronTablesMatPath = fullfile(outFolder, [prefix, '_', 'normalizedChevronTables.mat']);
measureTimeTablesMatPath = fullfile(outFolder, [prefix, '_', 'measureTimeTables.mat']);
popAvgTablesMatPath = fullfile(outFolder, [prefix, '_', 'popAvgTables.mat']);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
