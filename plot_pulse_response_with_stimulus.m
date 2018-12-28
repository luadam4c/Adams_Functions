function [fig, subplots, plots] = ...
                plot_pulse_response_with_stimulus (tVec, respVec, stimVec, varargin)
%% Plots a pulse response with its stimulus
% Usage: [fig, subplots, plots] = ...
%               plot_pulse_response_with_stimulus (tVec, respVec, stimVec, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       fig         - handle to figure
%                   specified as a figure object handle
%       subplots    - handles to subplots
%                   specified as a axes object handle array
%       plots       - handles to plots
%                   specified as a graphics object handle array
% Arguments:
%       tVec        - TODO: Description of reqarg1
%                   must be a numeric vector
%       respVec     - TODO: Description of reqarg1
%                   must be a numeric vector
%       stimVec     - TODO: Description of reqarg1
%                   must be a numeric vector
%       varargin    - 'OutFolder': the name of the directory for plots
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
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'SameAsPulse': whether always the same as 
%                                       the current pulse endpoints
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'MeanValueWindowMs': window in ms for 
%                                           calculating mean values
%                   must be a positive scalar
%                   default == 0.5 ms
%                   - 'MinPeakDelayMs': minimum peak delay (ms)
%                               after the end of the pulse
%                   must be a positive scalar
%                   default == 0 ms
%                   - 'ResponseParams': a table containing the parsed parameters
%                   must be a table containing at least the fields:
%                       baseValue
%                       steadyValue
%                       idxBaseStart
%                       idxBaseEnd
%                       idxSteadyStart
%                       idxSteadyEnd
%                   default == returned by parse_pulse_response.m
%                   - 'ParamsFile': a spreadsheet containing parsed parameters
%                   must be empty or a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       cd/annotation_in_plot.m
%       cd/compute_sampling_interval.m
%       cd/compute_axis_limits.m
%       cd/create_error_for_nargin.m
%       cd/isfigtype.m
%       cd/plot_horizontal_line.m
%       cd/plot_vertical_line.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/compute_and_plot_average_response.m
%       cd/compute_and_plot_all_responses.m

% File History:
% 2018-12-15 Moved from compute_and_plot_evoked_LFP.m
% 2018-12-17 Now uses compute_xlimits.m and compute_ylimits.m 
% 2018-12-17 Now plots dashed lines for baseValue and minPeakDelay
% 2018-12-17 Now uses unmatched varargin parts as parameters for plot()
% 2018-12-17 Now uses halfPeakValue and draws a double arrow for peak delay
% 2018-12-17 Now applies parse_pulse_response.m for defaults
%               if responseParams is not provided
% 2018-12-18 Now uses annotation_in_plot.m, plot_horizontal_line, 
%               plot_vertical_line.m
% 2018-12-18 Now can optionally read values from a 'ParamsFile'
% 2018-12-28 Now plots peak half width, peak 33% decay by double-exp fit
%               and peak 90% decay
% 

%% Hard-coded parameters
colorAnnotations = 'r';
colorLines = 'g';
labelsDefaultDefault = {'Response', 'Stimulus'};

%% Default values for optional arguments
labelsDefault = {};             % set later
outFolderDefault = pwd;         % save in present working directory by default
fileBaseDefault = 'Unnamed';
fileSuffixDefault = '_pulse_response';
responseNameDefault = 'Pulse Response';
saveFlagDefault = true;         % save the plot by default
figTypesDefault = 'png';        % default figure type(s) for saving
sameAsPulseDefault = true;      % use pulse endpoints by default
meanValueWindowMsDefault = 0.5; % calculating mean values over 0.5 ms by default
minPeakDelayMsDefault = 0;      % no minimum peak delay by default
responseParamsDefault = [];     % set later
paramsFileDefault = '';         % no parameters file by default

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
addParameter(iP, 'SameAsPulse', sameAsPulseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MeanValueWindowMs', meanValueWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addParameter(iP, 'ResponseParams', responseParamsDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'ParamsFile', paramsFileDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, tVec, respVec, stimVec, varargin{:});
labels = iP.Results.Labels;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
fileSuffix = iP.Results.FileSuffix;
responseName = iP.Results.ResponseName;
saveFlag = iP.Results.SaveFlag;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
sameAsPulse = iP.Results.SameAsPulse;
meanValueWindowMs = iP.Results.MeanValueWindowMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
responseParams = iP.Results.ResponseParams;
paramsFile = iP.Results.ParamsFile;

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

%% Preparation
% Parse responseParams if not provided
if isempty(responseParams)
    if ~isempty(paramsFile)
        responseParams = readtable(paramsFile);
    else
        % Compute the sampling interval
        siMs = compute_sampling_interval(tVec);

        % Parse the pulse response
        responseParams = parse_pulse_response(responseVecs, siMs, ...
                                    'PulseVector', pulseVector, ...
                                    'SameAsPulse', sameAsPulse, ...
                                    'MeanValueWindowMs', meanValueWindowMs, ...
                                    'MinPeakDelayMs', minPeakDelayMs);
    end
end

% Extract from responseParams
%   Note: must be consistent with parse_pulse_response.m
%   Note: this function assumes responseParams has only one row
minPeakDelayMs = responseParams.minPeakDelayMs;
idxBeforePulseEnd = responseParams.idxBeforePulseEnd;
baseValue = responseParams.baseValue;
minValueAfterMinDelay = responseParams.minValueAfterMinDelay;
maxValueAfterMinDelay = responseParams.maxValueAfterMinDelay;
idxPeak = responseParams.idxPeak;
indHalfWidthEnds = responseParams.indHalfWidthEnds;
idxPeakTimeConstant = responseParams.idxPeakTimeConstant;
idxPeak90Decay = responseParams.idxPeak90Decay;
peakValue = responseParams.peakValue;
peakAmplitude = responseParams.peakAmplitude;
halfPeakValue = responseParams.halfPeakValue;
peakTimeConstantValue = responseParams.peakTimeConstantValue;
peak90DecayValue = responseParams.peak90DecayValue;
peakDelayMs = responseParams.peakDelayMs;

% Extract from cell arrays
indHalfWidthEnds = indHalfWidthEnds{1};

% Decide on labels
if isempty(labels)
    % TODO: is_column_name.m and is_row_name.m for tables
    if ismember('labels', responseParams.Properties.VariableNames)
        labels = responseParams{:, 'labels'};

        % Make sure labels is a cell array of two elements
        if numel(labels) == 1
            labels = labels{1};
        end
    else
        labels = labelsDefaultDefault;
    end
end

% Check if needed output directories exist
check_dir(outFolder);

% Compute the y axis limits based on
%   the baseline value and the extrema after the minimum delay
[yLimitsResp, yRangeResp] = ...
    compute_axis_limits([minValueAfterMinDelay, maxValueAfterMinDelay, ...
                            baseValue], 'y');

% Compute the x axis limits
[xLimits, xRange] = compute_axis_limits(tVec, 'x');

% Compute times from indices
timePulseEnd = ...
    argfun(@(x) tVec(x), idxBeforePulseEnd);

% Compute times from delays
minPeakTime = timePulseEnd + minPeakDelayMs;

% If a peak exists, compute peak-related values
if ~isnan(idxPeak)
    % Compute the times relative to the minimum x value
    [timePeakRel, timePulseEndRel, timeHalfWidthEndsRel, ...
        timePeakTimeConstantRel, timePeak90DecayRel] = ...
        argfun(@(x) (tVec(x) - xLimits(1)) / xRange, ...
                idxPeak, idxBeforePulseEnd, indHalfWidthEnds, ...
                idxPeakTimeConstant, idxPeak90Decay);

    % Compute the values relative to the minimum y value
    [baseValueRel, peakValueRel, halfPeakValueRel, ...
        peakTimeConstantValueRel, peak90DecayValueRel] = ...
        argfun(@(x) (x - yLimitsResp(1)) / yRangeResp, ...
                baseValue, peakValue, halfPeakValue, ...
                peakTimeConstantValue, peak90DecayValue);

    % Compute the peak amplitude double arrow x and y values 
    %   in normalized units relative to the axes
    peakAmpXValues = timePeakRel * ones(1, 2);
    peakAmpYValues = sort([baseValueRel, peakValueRel]);

    % Compute the peak delay double arrow x and y values 
    %   in normalized units relative to the axes
    peakDelayXValues = [timePulseEndRel, timePeakRel];
    peakDelayYValues = baseValueRel * ones(1, 2);

    % Compute the peak half width double arrow x and y values 
    %   in normalized units relative to the axes
    peakHalfWidthXValues = timeHalfWidthEndsRel;
    peakHalfWidthYValues = halfPeakValueRel * ones(1, 2);

    % Compute the peak time constant double arrow x and y values 
    %   in normalized units relative to the axes
    peakTimeConstantXValues = [timePeakRel, timePeakTimeConstantRel];
    peakTimeConstantYValues = peakTimeConstantValue * ones(1, 2);

    % Compute the peak 90% decay double arrow x and y values 
    %   in normalized units relative to the axes
    peak90DecayXValues = [timePeakRel, timePeak90DecayRel];
    peak90DecayYValues = peak90DecayValue * ones(1, 2);

    % Create a labels for the peak amplitude and delay
    peakAmpLabel = ['peak amp = ', num2str(peakAmplitude)];
    peakDelayLabel = ['peak delay = ', num2str(peakDelayMs)];
end

%% Do the job
% Open and clear figure
if saveFlag
    fig = figure('Visible', 'off');
    figName = fullfile(outFolder, [fileBase, fileSuffix]);
    clf(fig);
else
    figure;
end

%% Generate a subplot for the pulse response
%   Annotations:
%       red double arrow for peak amplitude
subplots(1) = subplot(3, 1, 1:2); hold on;

% Plot the pulse response
plots(1) = plot(tVec, respVec, otherArguments);

% Update x-axis and y-axis limits
xlim(xLimits);
ylim(yLimitsResp);

% Plot a dashed horizontal line for baseValue
plots(2) = plot_horizontal_line(baseValue, 'XLimits', xLimits, ...
        'LineStyle', '--', 'Color', colorLines);

% Plot a dashed vertical line for timePulseEnd
plots(3) = plot_vertical_line(timePulseEnd, 'YLimits', yLimitsResp, ...
        'LineStyle', '--', 'Color', colorLines);

% Plot a dashed vertical line for minPeakTime
plots(4) = plot_vertical_line(minPeakTime, 'YLimits', yLimitsResp, ...
        'LineStyle', '--', 'Color', colorLines);

if ~isnan(idxPeak)
    % Draw a double arrow spanning the peak amplitude
    plots(5) = annotation_in_plot('doublearrow', peakAmpXValues, ...
                                peakAmpYValues, 'Color', colorAnnotations);

    % Draw a double arrow spanning the peak delay
    plots(6) = annotation_in_plot('doublearrow', peakDelayXValues, ...
                                peakDelayYValues, 'Color', colorAnnotations);

    % Draw a double arrow spanning the peak half width
    plots(7) = annotation_in_plot('doublearrow', peakHalfWidthXValues, ...
                            peakHalfWidthYValues, 'Color', colorAnnotations);

    % Draw a double arrow spanning the peak time constant
    plots(8) = annotation_in_plot('doublearrow', peakTimeConstantXValues, ...
                            peakTimeConstantYValues, 'Color', colorAnnotations);

    % Draw a double arrow spanning the peak 90% decay
    plots(9) = annotation_in_plot('doublearrow', peak90DecayXValues, ...
                            peak90DecayYValues, 'Color', colorAnnotations);

    % Show a text for the value of the peak amplitude above the peak
    text(tVec(idxPeak) - peakDelayMs * 1/4, ...
            peakValue + peakAmplitude * 1/16, peakAmpLabel);

    % Show a text for the value of the peak delay
    % TODO: use extract_elements(tVec, '1q', 'Endpoints', [idxBeforePulseEnd, idxPeak])
    text(tVec(round(idxBeforePulseEnd * 3/4 + idxPeak * 1/4)), ...
            baseValue - peakAmplitude * 1/16, peakDelayLabel);
end

% Generate a y-axis label
ylabel(labels{1});

% Generate a title for the pulse response
title([responseName, ' for ', fileBase], 'Interpreter', 'none');

%% Generate a subplot for the stimulation pulse
subplots(2) = subplot(3, 1, 3); hold on;

% Plot the stimulation pulse
plot(tVec, stimVec, otherArguments);

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
linkaxes([subplots(1), subplots(2)], 'x');

%% Save and close figure
if saveFlag
    save_all_figtypes(fig, figName, figTypes);
    close(fig)
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

minPeakDelayMsDefault = [];         % set later
idxBeforePulseEndDefault = [];      % set later
baseValueDefault = [];              % set later
minValueAfterMinDelayDefault = [];  % set later
maxValueAfterMinDelayDefault = [];  % set later
idxPeakDefault = [];                % set later
peakValueDefault = [];              % set later
halfPeakValueDefault = [];          % set later
peakAmplitudeDefault = [];          % set later
peakDelayMsDefault = [];            % set later
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'IdxBeforePulseEnd', idxBeforePulseEndDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
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
addParameter(iP, 'HalfPeakValue', halfPeakValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PeakAmplitude', peakAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PeakDelayMs', peakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
minPeakDelayMs = iP.Results.MinPeakDelayMs;
idxBeforePulseEnd = iP.Results.IdxBeforePulseEnd;
baseValue = iP.Results.BaseValue;
minValueAfterMinDelay = iP.Results.MinValueAfterMinDelay;
maxValueAfterMinDelay = iP.Results.MaxValueAfterMinDelay;
idxPeak = iP.Results.IdxPeak;
peakValue = iP.Results.PeakValue;
halfPeakValue = iP.Results.HalfPeakValue;
peakAmplitude = iP.Results.PeakAmplitude;
peakDelayMs = iP.Results.PeakDelayMs;

peakAmpLabel = ['peak amp = ', num2str(peakAmplitude), ' mV'];
peakDelayLabel = ['peak delay = ', num2str(peakDelayMs), ' ms'];

plots(2) = line(xLimits, baseValue * ones(size(xLimits)), ...
            'LineStyle', '--', 'Color', colorLines);

% Use plot_vertical_line.m
line(minPeakTime * ones(size(yLimitsResp)), yLimitsResp, ...
    'LineStyle', '--', 'Color', colorLines);

% Find the minimum and maximum values of interest,
%   including the baseline value and the trace after the minimum delay
minValueOfInterest = min([minValueAfterMinDelay, baseValue]);
maxValueOfInterest = max([maxValueAfterMinDelay, baseValue]);

%       cd/compute_xlimits.m
%       cd/compute_ylimits.m

% Compute appropriate y axis limits from the range of values of interest
[yLimitsResp, yRangeResp] = ...
    compute_ylimits(minValueOfInterest, maxValueOfInterest, 'Coverage', 80);

[xLimits, xRange] = compute_xlimits(tVec, 'Coverage', 100);

minPeakTime = xLimits(1) + minPeakDelayMs;

compute_axis_limits([maxValueAfterMinDelay, baseValue], 'y');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%