function m3ha_initial_slopes (varargin)
%% Computes and plots histograms of slopes right after current pulse start and end
% Usage: m3ha_initial_slopes (varargin)
% Arguments:
%       varargin    - 'HomeDir': the name of the directory 
%                                containing a /matfiles/ directory
%                   must be a directory
%                   default == '/media/adamX/m3ha/data_dclamp/take4'
%                   - 'RangeNSamples': the range of the number of samples 
%                       for taking the average slope
%                   must be a positive integer vector
%                   default == [2, 11]
%                   - 'HistFolder': the name of the directory where the 
%                                   histogram will be stored
%                   must be a directory
%                   default == fullfile(homeDir, 'histograms_initial_slopes')
%                   - 'CprFolder': the name of the directory where the sample 
%                                   cpr traces will be stored
%                   must be a directory
%                   default == fullfile(homeDir, 'cprtraces')
%                   - 'NCprToPlot': the number of current pulse responses to plot
%                                   at slope extremes
%                   must be a positive integer scalar
%                   default == 20
%                   - 'NSamplesForPlot': the number of samples to average when 
%                                   plotting CPR
%                   must be a positive integer scalar
%                   default == 10
%                   - 'ThresMethod': threshold method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'kernelVoid'    - use a kernel distribution
%                       'threeStdMainComponent' 
%                                       - use the center Gaussian distribution
%                   default == 'threeStdMainComponent'
%                   - 'UseCurrentFlag': whether to use the current trace
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       /media/adamX/m3ha/data_dclamp/take4/dclampdatalog_take4.mat
%       cd/check_dir.m
%       cd/correct_unbalanced_bridge.m
%       cd/dlmwrite_with_header.m
%       cd/find_initial_slopes.m
%       cd/fitdist_initial_slopes.m
%       cd/plot_histogram.m
%       cd/plot_pdf.m

% File History:
% 2018-05-24 BT - Created
% 2018-05-25 AL - Made many modifications
% 2018-05-29 BT - Added input parser
% 2018-05-29 AL - Modified initial slopes
% 2018-05-29 AL - Added use of log_arraytext
% 2018-06-04 BT - Fixed current pulse response to be vvecCpr;
%                   added current pulse response plotting
% 2018-06-05 AL - Finished Brian's modifications
% 2018-06-05 BT - Made 'NCprToPlot' an optional parameter
% 2018-06-05 AL - Now uses first and last sample points only for average slope
% 2018-06-05 BT - Made 'NSamplesForPlot' an optional parameter and updated cpr plotting
% 2018-06-06 AL - Changed abs(mean()) -> mean(abs())
% 2018-06-06 AL - Fixed usage of nSamplesForPlot
% 2018-06-06 AL - Fixed slope calculation
% 2018-06-06 AL - Now finds tvecCpr and vvecCpr only once
% 2018-06-11 AL - Now reverses the sign of the startSlope instead of using abs()
% 2018-06-11 AL - Now tries to account for systematic errors by trying
%                   pairs of points shifted to the right too
% 2018-06-15 AL - Changed histFolder to histograms_initial_slopes
% 2018-06-20 AL - Renamed m3ha_initial_slopes.m -> find_initial_slopes.m
% 2018-06-20 AL - Now saves variables in .mat files
%                   and changed the default nSamplesForPlot 10 -> 2
% 2018-06-21 AL - Now uses check_dir.m
% 2018-07-18 BT - Plots fitting with Gaussians
% 2018-07-18 AL - Revert back to using just the first pair of points
% 2018-07-18 AL - Order by absolute value of jump
% 2018-07-22 AL - Added peak3s
% 2018-07-22 AL - Now plots CPR traces near peak1 too
% 2018-08-11 AL - Added UseCurrentFlag and set the default to not use it
% 2018-08-12 AL - Set default of UseCurrentFlag to true
% 2018-08-12 AL - Now saves final output files for different thresMethod separately
% 2018-08-12 AL - Set default for thresMethod to threeStdMainComponent
% 2018-08-12 AL - Now fileNames is always a column vector
% 2018-08-12 AL - Now saves BestModels
% 2018-08-13 AL - Now saves all corrected traces too
% 2018-09-11 AL - Renamed find_initial_slopes.m -> m3ha_initial_slopes.m
% 2018-09-12 AL - Changed maxNumComponents from 4 -> 3

%% Constants
PA_PER_NA = 1000;                   % 1000 pA per nA

%% Hard-coded parameters (must be consistent with dataDclampExtractor.m)
dataType = 'mat';
cpStartWin = [95, 105];             % window of current pulse start (ms)
cpEndWin = [105, 115];              % window of current pulse end (ms)
outlierMethod = 'fiveStds';         % method for determining outliers
bandwidth2StdRatio = 1/5;           % ratio of kernel bandwidth to 
                                    %   standard deviation
% prec = 10^-4;   % Precision (must match fit_gaussians_and_refine_threshold.m)
validMethods = {'kernelVoid', 'threeStdMainComponent'};
maxNumComponents = 3; %4               % maximum number of components
                                    %   if multiple Gaussians are used to fit

%% Default values for optional arguments
rangeNSamplesDefault = [2, 11];
% homeDirDefault = '/media/adamX/m3ha/data_dclamp/take4'; % HomeDir on Linux
homeDirDefault = '/tmp/data/m3ha/data_dclamp/take4'; % HomeDir on Linux
% homeDirDefault = '/m3ha/data_dclamp/take4'; % HomeDir on Windows
cprFolderDefault = [];
histFolderDefault = [];
nCprToPlotDefault = 20;
nSamplesForPlotDefault = 2;
thresMethodDefault = 'threeStdMainComponent';
useCurrentFlagDefault = true;      % use the current trace by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Paths to functions on Windows
% addpath(fullfile('\\128.143.17.146\MatlabFishFish', 'Adams_Functions'))
% addpath(fullfile('\\128.143.17.146\MatlabFishFish', 'Downloaded_Functions'))
% addpath(fullfile('\\128.143.17.146\MatlabFishFish', 'Brians_Functions'))

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'HomeDir', homeDirDefault, @isdir); 
addParameter(iP, 'RangeNSamples', rangeNSamplesDefault, ...       
    @(x) validateattributes(x, {'numeric'}, ...
        {'vector', 'positive', 'integer', 'increasing', 'numel', 2}));
addParameter(iP, 'CprFolder', cprFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'HistFolder', histFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'NCprToPlot', nCprToPlotDefault, ...       
    @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'positive', 'integer'}));
addParameter(iP, 'NSamplesForPlot', nSamplesForPlotDefault, ...       
    @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'positive', 'integer'}));
addParameter(iP, 'ThresMethod', thresMethodDefault, ...
    @(x) any(validatestring(x, validMethods)));
addParameter(iP, 'UseCurrentFlag', useCurrentFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
rangeNSamples = iP.Results.RangeNSamples;
homeDir = iP.Results.HomeDir;
cprFolder = iP.Results.CprFolder;
histFolder = iP.Results.HistFolder;
nCprToPlot = iP.Results.NCprToPlot;
nSamplesForPlot = iP.Results.NSamplesForPlot;
thresMethod = validatestring(iP.Results.ThresMethod, validMethods);
useCurrentFlag = iP.Results.UseCurrentFlag;

% Set dependent argument defaults
if isempty(cprFolder)
    cprFolder = fullfile(homeDir, ...
                        ['cprtraces_nSamplesForPlot_', ...
                            num2str(nSamplesForPlot), '_', thresMethod]);
end
cprCorrectedFolder = fullfile(homeDir, ...
                    ['cprCorrectedtraces_nSamplesForPlot_', ...
                        num2str(nSamplesForPlot), '_', thresMethod]);
if isempty(histFolder)
    histFolder = fullfile(homeDir, ['histograms_initial_slopes_', thresMethod]);
end
histNewFolder = fullfile(homeDir, ['histograms_initial_slopes_', ...
                                    thresMethod, '_new']);

% Check relationships between arguments
if nSamplesForPlot < rangeNSamples(1) || ...
    nSamplesForPlot > rangeNSamples(2)
    error('nSamplesForPlot must be within rangeNSamples!');
end

%% Check if needed output directories exist
check_dir({histFolder, histNewFolder, cprFolder, cprCorrectedFolder});

%% Do the job
% Get the full path for the .mat files folder
matFolder = fullfile(homeDir, 'matfiles');

% Find all the files in the matfiles directory
dataFiles = dir(fullfile(matFolder, ['*.', dataType]));

% Find total number of sweeps to import
nSwps = numel(dataFiles);

% Extract the file names as a column cell array
dataFileNames = {dataFiles.name};
dataFileNames = dataFileNames(:);

% Convert to full file names
dataFullFileNames = ...
    cellfun(@(x) fullfile(matFolder, x), dataFileNames, 'UniformOutput', false);

% Get all possible nSamples to use
allNSamples = rangeNSamples(1):rangeNSamples(2);

% Count the number of possible nSamples to use
nNSamples = length(allNSamples);

%% Extract needed vectors
% Start timer
tic;

% Initialize cell array to store time and voltage traces
tvecCprAll = cell(nSwps, 1);
vvecCprAll = cell(nSwps, 1);
ivecCprAll = cell(nSwps, 1);

% Extract vectors
parfor iSwp = 1:nSwps
    % Get the current file name
    fileName = dataFullFileNames{iSwp};

    % Read in the matfile
    m = matfile(fileName);
    % fprintf('Using trace %s\n', fileName);

    % Use the time vector from original data trace (not resampled)
    tvecOrig = m.d_orig(:, 1);

    % Find the indices corresponding to the approximate current pulse window
    approxCpWinBegin = find(tvecOrig >= cpStartWin(1), 1);
    approxCpWinEnd = find(tvecOrig <= cpEndWin(2), 1, 'last');
    indApproxCpWin = approxCpWinBegin:approxCpWinEnd;

    % Restrict the time vector to current pulse window
    tvecCpr = tvecOrig(indApproxCpWin);

    % Restrict the current vector (in pA) to current pulse window
    ivecCpr = m.d_orig(indApproxCpWin, 3);

    % Convert the current vector to nA
    ivecCpr = ivecCpr / PA_PER_NA;

    % Restrict the voltage vector (in mV) to current pulse window
    vvecCpr = m.d_orig(indApproxCpWin, 4);

    % Store the time and voltage traces for plotting
    tvecCprAll{iSwp} = tvecCpr;
    vvecCprAll{iSwp} = vvecCpr;
    ivecCprAll{iSwp} = ivecCpr;
end

% End timer
toc;

%% Find and log all initial slopes
% Start timer
tic;

% Create file names
csvFileName = fullfile(homeDir, 'initial_slopes_avgSlopes.csv');
matFile1Name = fullfile(homeDir, 'initial_slopes_avgSlopes.mat');
matFile2Name = fullfile(homeDir, 'initial_slopes_for_plotting.mat');

% Find and log all initial slopes
[startSlopes, endSlopes, avgSlopes, isUnbalancedAll, indUsedForPlot] = ...
    find_and_log_initial_slopes(tvecCprAll, ivecCprAll, ...
                                vvecCprAll, allNSamples, nNSamples, ...
                                csvFileName, matFile1Name, matFile2Name, ...
                                nSamplesForPlot, useCurrentFlag);

% End timer
toc;

%% Fit the initial slope distributions
% Start timer
tic;

% Load average slopes from a matfile
load(fullfile(homeDir, 'initial_slopes_avgSlopes.mat'), ...
                        'allNSamples', 'avgSlopes');

% Fit all distributions
[bestModels, pdfModels, peak1s, peak2s, peak3s, threshold1s, threshold2s] = ...
    fitdist_initial_slopes(allNSamples, avgSlopes, ...
                            'ThresMethod', thresMethod, ...
                            'OutlierMethod', outlierMethod, ...
                            'Bandwidth2StdRatio', bandwidth2StdRatio, ...
                            'MaxNumComponents', maxNumComponents, ...
                            'OutFolder', histFolder);

% End timer
toc;

%% Correct traces
%   Note: This assumes nSamples = 2
vvecCprCorrectedAll = cell(nSwps, 1);
parfor iSwp = 1:nSwps
    vvecCprCorrectedAll{iSwp} = ...
        correct_unbalanced_bridge(vvecCprAll{iSwp}, ivecCprAll{iSwp}, ...
                                  'UseCurrentFlag', useCurrentFlag);
end

%% Sort current pulse responses by average slopes and 
%   determine the new vvecCpr to be used

% Determine the index of nSamplesForPlot in allNSamples
idxNSamplesForPlot = find(allNSamples == nSamplesForPlot);

% Sort the average slopes for the nSamples value to plot in descending order
[avgSlopesSorted, origInd] = sort(avgSlopes(:, idxNSamplesForPlot), 'descend');

% Reorder everything else the same way
filenamesSorted = dataFullFileNames(origInd);
tvecCprAllSorted = tvecCprAll(origInd);
ivecCprAllSorted = ivecCprAll(origInd);
vvecCprAllSorted = vvecCprAll(origInd);
vvecCprCorrectedAllSorted = vvecCprCorrectedAll(origInd);
indUsedForPlotSorted = indUsedForPlot(origInd, :);
startSlopesSorted = startSlopes(origInd, idxNSamplesForPlot);
endSlopesSorted = endSlopes(origInd, idxNSamplesForPlot);
isUnbalancedAllSorted = isUnbalancedAll(origInd, idxNSamplesForPlot);

% Find the thresholds for the nSamples value to plot
threshold1ForPlot = threshold1s(idxNSamplesForPlot);
peak1ForPlot = peak1s(idxNSamplesForPlot);
threshold2ForPlot = threshold2s(idxNSamplesForPlot);

% Find the iSort for the current pulse response near the thresholds and peak1
iThreshold1Balanced = find(avgSlopesSorted < threshold1ForPlot, 1, 'first');
iPeak1Balanced = find(avgSlopesSorted < peak1ForPlot, 1, 'first');
iThreshold2Balanced = find(avgSlopesSorted > threshold2ForPlot, 1, 'last');

% Determine vvecCprNew: Use the corrected trace for those traces
%   with average slopes outside the thresholds; use the original trace otherwise
vvecCprNewSorted = vvecCprCorrectedAllSorted;
vvecCprNewSorted(iThreshold2Balanced:iThreshold1Balanced) = ...
    vvecCprAllSorted(iThreshold2Balanced:iThreshold1Balanced);

%% Find and log all initial slopes again
% Start timer
tic;

% Create file names
csvFileName = fullfile(homeDir, 'initial_slopes_avgSlopes_new.csv');
matFile1Name = fullfile(homeDir, 'initial_slopes_avgSlopes_new.mat');
matFile2Name = fullfile(homeDir, 'initial_slopes_for_plotting_new.mat');

% Find and log all initial slopes
[startSlopesNew, endSlopesNew, avgSlopesNew, ...
    isUnbalancedAllNew, indUsedForPlotNew] = ...
    find_and_log_initial_slopes(tvecCprAllSorted, ivecCprAllSorted, ...
                                vvecCprNewSorted, allNSamples, nNSamples, ...
                                csvFileName, matFile1Name, matFile2Name, ...
                                nSamplesForPlot, useCurrentFlag);

% End timer
toc;

%% Fit the initial slope distributions again
% Start timer
tic;

% Fit all distributions
[bestModelsNew, pdfModelsNew, peak1sNew, peak2sNew, peak3sNew, ...
    threshold1sNew, threshold2sNew] = ...
    fitdist_initial_slopes(allNSamples, avgSlopesNew, ...
                            'ThresMethod', thresMethod, ...
                            'OutlierMethod', outlierMethod, ...
                            'Bandwidth2StdRatio', bandwidth2StdRatio, ...
                            'MaxNumComponents', maxNumComponents, ...
                            'OutFolder', histNewFolder);

% End timer
toc;

%% Store all variables in a matfile
% Create mat file name
matFileName = fullfile(homeDir, ['initial_slopes_nSamplesForPlot_', ...
                                    num2str(nSamplesForPlot), ...
                                    '_', thresMethod, '.mat']);

% Save sorted variables to a matfile
save(matFileName, ...
    'bestModels', 'pdfModels', 'peak1s', 'peak2s', 'peak3s', ...
    'threshold1s', 'threshold2s', ...
    'allNSamples', 'nSamplesForPlot', ...
    'origInd', 'threshold1ForPlot', 'threshold2ForPlot', ...
    'iThreshold1Balanced', 'iThreshold2Balanced', ...
    'filenamesSorted', 'indUsedForPlotSorted', ...
    'tvecCprAllSorted', 'vvecCprAllSorted', 'vvecCprNewSorted', ...
    'startSlopesSorted', 'endSlopesSorted', 'avgSlopesSorted', ...
    'isUnbalancedAllSorted', 'vvecCprCorrectedAll', ...
    'startSlopesNew', 'endSlopesNew', 'avgSlopesNew', ...
    'bestModelsNew', 'pdfModelsNew', 'peak1sNew', 'peak2sNew', 'peak3sNew', ...
    'threshold1sNew', 'threshold2sNew', ...
    '-v7.3');
%    'allNSamples', 'csvHeader', 'nSamplesForPlot', ...

%% Plot histograms for each nSample
% Start timer
tic;

combineOutliersFlag = true;
plot_initial_slope_histograms(avgSlopes, pdfModels, ...
                                peak1s, peak2s, peak3s, ...
                                threshold1s, threshold2s, ...
                                histFolder, nNSamples, allNSamples, ...
                                outlierMethod, combineOutliersFlag);
combineOutliersFlag = false;
plot_initial_slope_histograms(avgSlopes, pdfModels, ...
                                peak1s, peak2s, peak3s, ...
                                threshold1s, threshold2s, ...
                                histFolder, nNSamples, allNSamples, ...
                                outlierMethod, combineOutliersFlag);
combineOutliersFlag = true;
plot_initial_slope_histograms(avgSlopesNew, pdfModelsNew, ...
                                peak1sNew, peak2sNew, peak3sNew, ...
                                threshold1sNew, threshold2sNew, ...
                                histNewFolder, nNSamples, allNSamples, ...
                                outlierMethod, combineOutliersFlag);
combineOutliersFlag = false;
plot_initial_slope_histograms(avgSlopesNew, pdfModelsNew, ...
                                peak1sNew, peak2sNew, peak3sNew, ...
                                threshold1sNew, threshold2sNew, ...
                                histNewFolder, nNSamples, allNSamples, ...
                                outlierMethod, combineOutliersFlag);

% End timer
toc;

%% Plot sample current pulse responses with largest average slopes
% Start timer
tic;

% Load variables used for plotting
load(fullfile(homeDir, 'initial_slopes_for_plotting.mat'), ...
        'tvecCprAll', 'vvecCprAll', 'indUsedForPlot');

% Find the iSort of interest for plotting
iSortMostPositive = 1:nCprToPlot;
if isempty(iThreshold1Balanced)
    iSortNearThreshold1 = [];
else
    iSortNearThreshold1 = (iThreshold1Balanced - floor(nCprToPlot/2) + 1): ...
                            (iThreshold1Balanced + ceil(nCprToPlot/2));
end
iSortNearPeak1 = (iPeak1Balanced - floor(nCprToPlot/2) + 1): ...
                        (iPeak1Balanced + ceil(nCprToPlot/2));
if isempty(iThreshold2Balanced)
    iSortNearThreshold2 = [];
else
    iSortNearThreshold2 = (iThreshold2Balanced - floor(nCprToPlot/2) + 1): ...
                            (iThreshold2Balanced + ceil(nCprToPlot/2));
end
iSortMostNegative = (nSwps-nCprToPlot+1):nSwps;
iSortOfInterest = ...
    union(iSortMostPositive, ...
        union(iSortNearThreshold1, ...
            union(iSortNearPeak1, ...
                union(iSortNearThreshold2, iSortMostNegative))));

% Remove any indices out of range
iSortOfInterest = ...
    iSortOfInterest(iSortOfInterest >= 1 & iSortOfInterest <= nSwps);

% Plot the current pulse responses of interest
parfor j = 1:length(iSortOfInterest)
    % Get the current iSort
    iSort = iSortOfInterest(j);

    % Plot the corresponding current pulse response
    plot_cpr(iSort, iThreshold1Balanced, iThreshold2Balanced, ...
             cprFolder, ...
             tvecCprAllSorted, vvecCprAllSorted, [], ...
             indUsedForPlotSorted, filenamesSorted, ...
             startSlopesSorted, endSlopesSorted, avgSlopesSorted, ...
             isUnbalancedAllSorted);
    plot_cpr(iSort, iThreshold1Balanced, iThreshold2Balanced, ...
             cprCorrectedFolder, ...
             tvecCprAllSorted, vvecCprAllSorted, vvecCprCorrectedAllSorted, ...
             indUsedForPlotSorted, filenamesSorted, ...
             startSlopesSorted, endSlopesSorted, avgSlopesSorted, ...
             isUnbalancedAllSorted); 
end

% Plot all current pulse responses if nSamplesForPlot == 2
if nSamplesForPlot == 2 && strcmp(thresMethod, 'threeStdMainComponent')
    % Construct a separate outfolder for all current pulse response traces
    cprAllFolder = [cprFolder, '_all'];
    cprCorrectedAllFolder = [cprCorrectedFolder, '_all'];
    check_dir({cprAllFolder, cprCorrectedAllFolder});

    parfor iSort = 1:nSwps
        plot_cpr(iSort, iThreshold1Balanced, iThreshold2Balanced, ...
                 cprAllFolder, ...
                 tvecCprAllSorted, vvecCprAllSorted, [], ...
                 indUsedForPlotSorted, filenamesSorted, ...
                 startSlopesSorted, endSlopesSorted, avgSlopesSorted, ...
                 isUnbalancedAllSorted);
        plot_cpr(iSort, iThreshold1Balanced, iThreshold2Balanced, ...
                 cprCorrectedAllFolder, ...
                 tvecCprAllSorted, vvecCprAllSorted, vvecCprCorrectedAllSorted, ...
                 indUsedForPlotSorted, filenamesSorted, ...
                 startSlopesSorted, endSlopesSorted, avgSlopesSorted, ...
                 isUnbalancedAllSorted);
    end
end

% End timer
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [startSlopes, endSlopes, avgSlopes, ...
            isUnbalancedAll, indUsedForPlot] = ...
    find_and_log_initial_slopes(tvecCprAll, ivecCprAll, ...
                                vvecCprAll, allNSamples, nNSamples, ...
                                csvFileName, matFile1Name, matFile2Name, ...
                                nSamplesForPlot, useCurrentFlag)

%% Find all initial slopes
[startSlopes, endSlopes, avgSlopes, isUnbalancedAll, indUsedForPlot] = ...
    find_initial_slopes(tvecCprAll, ivecCprAll, vvecCprAll, allNSamples, ...
                        'NSamplesForPlot', nSamplesForPlot, ...
                        'UseCurrentFlag', useCurrentFlag);

%% Log average slopes in a CSV file
% Construct header for csv file
csvHeader = cellfun(@(x) ['nSamples', num2str(x)], ...
                    num2cell(allNSamples), 'UniformOutput', false);

% Write csv file with header
dlmwrite_with_header(csvFileName, avgSlopes, 'ColumnHeader', csvHeader);

%% Save variables in matfiles
% Save average slopes in a matfile
save(matFile1Name, ...
               'allNSamples', 'nNSamples', ...
               'startSlopes', 'endSlopes', 'avgSlopes', ...
               'isUnbalancedAll', '-v7.3');


% Save variables used for plotting in a matfile
save(matFile2Name, ...
               'tvecCprAll', 'vvecCprAll', 'indUsedForPlot', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_initial_slope_histograms(avgSlopes, pdfModels, ...
                                    peak1s, peak2s, peak3s, ...
                                    threshold1s, threshold2s, ...
                                    histFolder, nNSamples, allNSamples, ...
                                    outlierMethod, combineOutliersFlag)

% Use histogram()
parfor iNSample = 1:nNSamples
    % Get the current nSamples
    nSamples = allNSamples(iNSample);

    % Get the current average slopes
    avgSlopesThis = avgSlopes(:, iNSample);

    % Remove outliers if any
    avgSlopesThisTrunc = ...
        remove_outliers(avgSlopesThis, 'OutlierMethod', outlierMethod);

    % Construct a table for the lines to plot
    Value = [peak1s(iNSample); peak2s(iNSample); peak3s(iNSample); ...
             threshold1s(iNSample); threshold2s(iNSample)];
    Color = {'Red'; 'SpringGreen'; 'DarkGreen'; 'DeepSkyBlue'; 'MediumBlue'};
    LineStyle = {'--'; '--'; '--'; '--'; '--'};
    Label = {['Peak 1 = ', num2str(peak1s(iNSample)), ' (V/s)']; ...
             ['Peak 2 = ', num2str(peak2s(iNSample)), ' (V/s)']; ...
             ['Peak 3 = ', num2str(peak3s(iNSample)), ' (V/s)']; ...
             ['Threshold 1 = ', num2str(threshold1s(iNSample)), ' (V/s)']; ...
             ['Threshold 2 = ', num2str(threshold2s(iNSample)), ' (V/s)']};
    linesToPlot = table(Value, Color, LineStyle, Label);

    % Plot the full histogram
    h = figure('Visible', 'off');
    clf(h);
    hold on
    if combineOutliersFlag
        plot_histogram(avgSlopesThis, 'OutlierMethod', outlierMethod);
        plot_pdf(avgSlopesThisTrunc, 'PDF', pdfModels{iNSample}, ...
                                     'LinesToPlot', linesToPlot);
        title(['Initial Slopes with rangeNSamples == ', ...
                num2str(nSamples), ' (zoomed)']);
        figname = fullfile(histFolder, ['initial_slope_nSamples_', ...
                num2str(nSamples), '_zoomed']);
    else
        histogram(avgSlopesThis);
        plot_pdf(avgSlopesThis, 'PDF', pdfModels{iNSample}, ...
                                'LinesToPlot', linesToPlot);
        title(['Initial Slopes for All Traces with rangeNSamples == ', ...
                num2str(nSamples)]);
        figname = fullfile(histFolder, ['initial_slope_nSamples_', ...
                num2str(nSamples)]);
    end
    xlabel('Slope (V/s)');
    ylabel('Trace Count');
    legend('location', 'northeast');
    saveas(h, figname, 'png');
    close(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_cpr(iSort, iThreshold1Balanced, iThreshold2Balanced, ...
                  cprFolder, tvecCprAllSorted, vvecCprAllSorted, ...
                  vvecCprCorrectedAllSorted, indUsedForPlotSorted, ...
                  filenamesSorted, startSlopesSorted, endSlopesSorted, ...
                  avgSlopesSorted, isUnbalancedAllSorted)

% Decide on the title color
if (~isempty(iThreshold2Balanced) && iSort < iThreshold2Balanced) || ...
    (~isempty(iThreshold1Balanced) && iSort > iThreshold1Balanced)
    titleColor = 'r';
else
    titleColor = 'k';
end

% Get the file name and file base
fileName = filenamesSorted{iSort};
[~, fileBase, ~] = fileparts(fileName);

% Get the time and voltage traces
tvecCpr = tvecCprAllSorted{iSort};
vvecCpr = vvecCprAllSorted{iSort};

% Get corrected voltage trace if using
if ~isempty(vvecCprCorrectedAllSorted) && strcmp(titleColor, 'r')
    vvecCprCorrected = vvecCprCorrectedAllSorted{iSort};
end

% Get the average slope
startSlope = startSlopesSorted(iSort);
endSlope = endSlopesSorted(iSort);
avgSlope = avgSlopesSorted(iSort);
indUsed = indUsedForPlotSorted(iSort, :);
isUnbalanced = isUnbalancedAllSorted(iSort);

% Create a figure
h = figure('Visible', 'off');
clf(h)
hold on

% Plot the original voltage trace in blue
plot(tvecCpr, vvecCpr, 'b.');

% Plot other traces and annotations
if ~isempty(vvecCprCorrectedAllSorted) && strcmp(titleColor, 'r')
    % Plot the corrected voltage trace in green
    plot(tvecCpr, vvecCprCorrected, 'g.');

    % Circle the points used for computing slopes in green or blue
    plot(tvecCpr(indUsed), vvecCprCorrected(indUsed), 'go');
    plot(tvecCpr(indUsed), vvecCpr(indUsed), 'bo');

    % Plot the slopes in green or blue
    line(tvecCpr(indUsed(1:2)), vvecCprCorrected(indUsed(1:2)), ...
            'Color', 'green');
    line(tvecCpr(indUsed(3:4)), vvecCprCorrected(indUsed(3:4)), ...
            'Color', 'green');
    line(tvecCpr(indUsed(1:2)), vvecCpr(indUsed(1:2)), 'Color', 'blue');
    line(tvecCpr(indUsed(3:4)), vvecCpr(indUsed(3:4)), 'Color', 'blue');
else
    % Circle the points used for computing slopes in red or blue
    if isUnbalanced
        plot(tvecCpr(indUsed), vvecCpr(indUsed), 'ro');
    else
        plot(tvecCpr(indUsed), vvecCpr(indUsed), 'bo');
    end

    % Plot the slopes in red
    line(tvecCpr(indUsed(1:2)), vvecCpr(indUsed(1:2)), 'Color', 'red');
    line(tvecCpr(indUsed(3:4)), vvecCpr(indUsed(3:4)), 'Color', 'red');
end

xlabel('Time (ms)');
ylabel('Voltage (mV)');
title(['Current Pulse Response for trace #', num2str(iSort), ...
        ' (', fileBase, ')'], 'Color', titleColor, 'Interpreter', 'none');
text(0.6, 0.95, ['Start slope (V/s) = ', num2str(startSlope)], ...
        'Units', 'normalized', 'Interpreter', 'none');
text(0.6, 0.9, ['End slope (V/s) = ', num2str(endSlope)], ...
        'Units', 'normalized', 'Interpreter', 'none');
text(0.6, 0.85, ['Average slope (V/s) = ', num2str(avgSlope)], ...
        'Units', 'normalized', 'Color', titleColor, 'Interpreter', 'none');
figname = fullfile(cprFolder, ['cpr_', num2str(iSort), '_', fileBase]);
saveas(h, figname, 'png');
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
