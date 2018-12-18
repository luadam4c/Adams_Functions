function h = plot_pulse_response_with_stimulus (tVec, respVec, stimVec, varargin)
%% Plots a pulse response with its stimulus
% Usage: h = plot_pulse_response_with_stimulus (tVec, respVec, stimVec, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - handle to figure
%                   specified as a figure handle
% Arguments:
%       tVec        - TODO: Description of reqarg1
%                   must be a numeric vector
%       respVec     - TODO: Description of reqarg1
%                   must be a numeric vector
%       stimVec     - TODO: Description of reqarg1
%                   must be a numeric vector
%       varargin    - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'OutFolder': the name of the directory for plots
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'SaveFlag': whether to save the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/compute_xlimits.m
%       cd/compute_ylimits.m
%       cd/create_error_for_nargin.m
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/compute_and_plot_average_response.m
%       cd/plot_protocols.m

% File History:
% 2018-12-15 Moved from compute_and_plot_evoked_LFP.m
% 2018-12-17 Now uses compute_xlimits.m and compute_ylimits.m 
% 2018-12-17 Now plots dashed lines for baseValue and minPeakDelay
% TODO: Deal with default behavior for parameters
% TODO: Read values from a sheetFile
% 

%% Hard-coded parameters
colorAnnotations = 'r';
colorLines = 'g';

%% Default values for optional arguments
minPeakDelayMsDefault = [];         % set later
baseValueDefault = [];              % set later
minValueAfterMinDelayDefault = [];  % set later
maxValueAfterMinDelayDefault = [];  % set later
idxPeakDefault = [];                % set later
peakValueDefault = [];              % set later
peakAmplitudeDefault = [];          % set later
peakDelayMsDefault = [];            % set later
labelsDefault = [];                 % set later
outFolderDefault = pwd;         % save in present working directory by default
fileBaseDefault = 'Unnamed';
fileSuffixDefault = '_pulse_response';
responseNameDefault = 'Pulse Response';
saveFlagDefault = true;         % save the plot by default
figTypesDefault = 'png';        % default figure type(s) for saving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'tVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'respVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'stimVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'BaseValue', baseValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MinValueAfterMinDelay', minValueAfterMinDelayDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MaxValueAfterMinDelay', maxValueAfterMinDelayDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'IdxPeak', idxPeakDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'PeakValue', peakValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PeakAmplitude', peakAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PeakDelayMs', peakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Labels', labelsDefault, ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileSuffix', fileSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ResponseName', responseNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, tVec, respVec, stimVec, varargin{:});
minPeakDelayMs = iP.Results.MinPeakDelayMs;
baseValue = iP.Results.BaseValue;
minValueAfterMinDelay = iP.Results.MinValueAfterMinDelay;
maxValueAfterMinDelay = iP.Results.MaxValueAfterMinDelay;
idxPeak = iP.Results.IdxPeak;
peakValue = iP.Results.PeakValue;
peakAmplitude = iP.Results.PeakAmplitude;
peakDelayMs = iP.Results.PeakDelayMs;
labels = iP.Results.Labels;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
fileSuffix = iP.Results.FileSuffix;
responseName = iP.Results.ResponseName;
saveFlag = iP.Results.SaveFlag;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Preparation
% Check if needed output directories exist
check_dir(outFolder);

% Compute appropriate y axis limits from the range of values 
%   after minimum peak delay
[yLimitsResp, yRangeResp] = ...
    compute_ylimits(minValueAfterMinDelay, maxValueAfterMinDelay);

% Compute the x axis limits
[xLimits, xRange] = compute_xlimits(tVec);

% Compute the time of the peak relative to the minimum x value
timePeakRel = (tVec(idxPeak) - xLimits(1)) / xRange;

% Compute the peak amplitude double array x and y values 
%   in normalized units relative to the subplot
peakAmpXValues = timePeakRel * ones(1, 2);
peakAmpYValues = ([baseValue, peakValue] - yLimitsResp(1)) / yRangeResp;

% Create a label for the peak amplitude
peakAmpLabel = ['peak amp = ', num2str(peakAmplitude), ' (mV)'];

%% Do the job
% Open and clear figure
if saveFlag
    h = figure('Visible', 'off');
    figName = fullfile(outFolder, [fileBase, fileSuffix]);
    clf(h);
else
    figure;
end

%% Generate a subplot for the pulse response
%   Annotations:
%       red double arrow for peak amplitude
ax1 = subplot(3, 1, 1:2); hold on;

% Plot the pulse response
plot(tVec, respVec);

% Update x-axis and y-axis limits
xlim(xLimits);
ylim(yLimitsResp);

% Compute the peak amplitude double array x and y values 
%   in normalized units relative to the figure
pos = get(gca, 'Position');
peakAmpXPositions = pos(1) + pos(3) * peakAmpXValues;
peakAmpYPositions = pos(2) + pos(4) * peakAmpYValues;

% Plot a dashed line for baseline
% Use plot_horizontal_line.m
line(tVec, xLimits, baseValue * ones(size(xLimits)), ...
    'LineStyle', '--', 'Color', colorLines);

% Plot a dashed line for minPeakDelay
% Use plot_vertical_line.m
line(tVec, minPeakDelayMs * ones(size(yLimitsResp)), yLimitsResp, ...
    'LineStyle', '--', 'Color', colorLines);

% Draw a doublearrow spanning the peak amplitude
annotation('doublearrow', peakAmpXPositions, peakAmpYPositions, ...
            'Color', colorAnnotations);

% Show a text for the value of the peak amplitude
text(tVec(idxPeak) + 1, mean([baseValue, peakValue]), peakAmpLabel);

% Generate a y-axis label
ylabel(labels{1});

% Generate a title for the pulse response
title([responseName, ' for ', fileBase], 'Interpreter', 'none');

%% Generate a subplot for the stimulation pulse
ax2 = subplot(3, 1, 3); hold on;

% Plot the stimulation pulse
plot(tVec, stimVec);

% Update x-axis limits
xlim(xLimits);    

% Generate a y-axis label
ylabel(labels{2});

% Generate an x-axis label
xlabel('Time (ms)');

% Generate a title for the stimulation pulse
title(['Stimulus for ', fileBase], 'Interpreter', 'none');

%% Finish up figure
% Link the axes on the subplots
linkaxes([ax1, ax2], 'x');

%% Save and close figure
if saveFlag
    save_all_figtypes(h, figName, figTypes);
    close(h)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

yLimitsResp = get(gca, 'YLim');

pos = get(gca, 'Position');
annotation('doublearrow', pos(1) + pos(3) * timePeakRel * ones(1, 2), ...
            pos(2) + pos(4) * ([baseValue, peakValue] - yLimitsResp(1)) / yRangeResp, ...
            'Color', colorAnnotations);

left = min(tVec);
right = max(tVec);
xLimits = [left, right];
% Compute the time range
xRange = right - left;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%