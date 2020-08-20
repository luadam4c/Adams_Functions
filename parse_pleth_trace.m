function varargout = parse_pleth_trace (vector, siSeconds, varargin)
%% Parses pleth traces
% Usage: [respParams, respData, handles] = ...
%                   parse_pleth_trace (vector, siSeconds, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/Data_for_test_analysis/parse_pleth_trace/20200811T175223_000.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%       channelValues = spike2Table.channelValues;
%       channelNames = spike2Table.channelNames;
%       plethVec = channelValues{strcmp(channelNames, 'Pleth 2')};
%       siSeconds = spike2Table{strcmp(channelNames, 'Pleth 2'), 'siSeconds'};
%       timeVecSec = create_time_vectors('TimeUnits', 's', 'SamplingIntervalSeconds', siSeconds, 'Vectors', plethVec);
%       [respParams, respData] = parse_pleth_trace(plethVec, siSeconds, 'TraceFileName', spike2MatPath, 'PlotPsdFlag', true, 'PlotRespFlag', true);
%
% Outputs:
%       respParams  - respiration measure parameters, including:
%                       traceFileName
%                       pathBase
%                       plethWindowSeconds
%                       resolutionSeconds
%                   specified as a scalar structure
%       respParams  - respiration measure data
%                       endPointsWindows
%                       timeInstantsSec
%                       respRates   - respiratory rates in rpm
%                       respAmps    - respiratory amplitudes in Volts
%                   specified as a scalar structure
%
% Arguments:
%       vector     - pleth trace
%                   TODO: Multiple traces
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vector
%       siSeconds   - sampling interval in seconds
%                   must be a positive vector
%       varargin    - 'TraceFileName': Name of the corresponding trace file(s)
%                   must be empty, a character vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - 'PathBase': file path base for output files
%                   must be a string scalar or a character vector
%                   default == match the trace path base
%                   - 'PlethWindowSeconds': pleth window for computing
%                                               breathing rate in seconds
%                   must be a numeric scalar
%                   default = 60 seconds
%                   - 'ResolutionSeconds': resolution in seconds
%                   must be a numeric scalar
%                   default = 60 seconds
%                   - 'PlotPsdFlag': whether to plot psd
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotRespFlag': whether to plot resp traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/argfun.m
%       cd/array_fun.m
%       cd/compute_running_windows.m
%       cd/convert_units.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_subplots.m
%       cd/create_time_vectors.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/freqfilter.m
%       cd/isemptycell.m
%       cd/match_row_count.m
%       cd/match_format_vector_sets.m
%       cd/parse_psd.m
%       cd/plot_traces.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/parse_spike2_mat.m

% File History:
% 2020-08-12 Modified from parse_laser_trace.m
% 2020-08-13 Now outputs a resp table
% 2020-08-14 Now restricts respRate to a range

%% Hard-coded parameters
% Decide on x and y limits to plot
xLimitsMinAll = {[20, 60]; [80, 120]; [140, 180]; ...
                [200, 240]; [260, 300]; [320, 360]};
yLimitsPleth = [-5, 5];
yLimitsRespRate = [];
yLimitsRespAmp = [];

% TODO: Make optional argument
nWindowsToPlot = 10;            % Number of windows to plot
respRateRangeRpm = [30, 150];   % Resp rate range in rpm

%% Default values for optional arguments
traceFileNameDefault = '';      % set later
plethWindowSecondsDefault = 60;
resolutionSecondsDefault = 60;
pathBaseDefault = '';
plotPsdFlagDefault = false;
plotRespFlagDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vector', ...                   % vector
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vector must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siSeconds', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TraceFileName', traceFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PathBase', pathBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlethWindowSeconds', plethWindowSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ResolutionSeconds', resolutionSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PlotPsdFlag', plotPsdFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotRespFlag', plotRespFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vector, siSeconds, varargin{:});
traceFileName = iP.Results.TraceFileName;
pathBase = iP.Results.PathBase;
plethWindowSeconds = iP.Results.PlethWindowSeconds;
resolutionSeconds = iP.Results.ResolutionSeconds;
plotPsdFlag = iP.Results.PlotPsdFlag;
plotRespFlag = iP.Results.PlotRespFlag;

%% Preparation
% Force vector to be column vector
vector = force_column_vector(vector);

% Count the number of samples
nSamples = count_samples(vector);

% Count the number of vector
nVectors = count_vectors(vector);

% Create time vector in seconds
timeVecSec = create_time_vectors(nSamples, 'TimeUnits', 's', ...
                    'SamplingIntervalSeconds', siSeconds, 'Vectors', vector);

% Match the row count
siSeconds = match_row_count(siSeconds, nVectors);

% Match the number of file names to the number of vector
[traceFileName, vector] = match_format_vector_sets(traceFileName, vector);

% Make sure pathBase is in agreement with nVectors
if isempty(pathBase)
    if isempty(traceFileName)
        pathBase = create_time_stamp;
    elseif isemptycell(traceFileName)
        pathBase = create_labels_from_numbers(1:nVectors, ...
                                'Prefix', strcat(create_time_stamp, '_'));
    else
        % Extract from the trace file name
        pathBase = extract_fileparts(traceFileName, 'pathbase');
    end
else
    if ischar(pathBase) && nVectors > 1 || ...
            iscell(pathBase) && numel(pathBase) ~= nVectors
        fprintf('Number of file bases must match the number of vector!\n');
        varargout{1} = table.empty;
        varargout{2} = table.empty;
        return
    end
end

%% Do the job
% Compute running windows
[endPointsWindows, timeInstantsSec] = ...
    compute_running_windows(timeVecSec, plethWindowSeconds, 'IsRegular', true, ...
                                'Resolution', resolutionSeconds);

% Count the number of windows
nWindows = size(endPointsWindows, 2);

% Create new file bases
pathBases = create_labels_from_numbers(1:nWindows, ...
                            'Prefix', strcat(pathBase, '_pleth_window'));

% Create window numbers
windowNumbers = (1:nWindows)';

% Compute base for window numbers
windowNumberBase = floor(nWindows/nWindowsToPlot);

% Decide whether to plot each window
plotFlagsEachWindow = plotPsdFlag & mod(windowNumbers, windowNumberBase) == 0;

% Compute respiratory rates and amplitudes for each running window
[respRates, respAmps] = ...
    array_fun(@(x, y) parse_pleth_trace_helper(vector, endPointsWindows(:, x), ...
                            siSeconds, pathBases{x}, y, respRateRangeRpm), ...
        	windowNumbers, plotFlagsEachWindow);

%% Plot results
if plotRespFlag
    % Convert times to minutes
    [timeVecsMin, timeInstantsMin] = ...
        argfun(@(x) convert_units(x, 's', 'min'), timeVecSec, timeInstantsSec);

    % Plot respiratory measure traces
    handles = cellfun(@(x) plot_resp_traces(timeVecsMin, vector, ...
                        timeInstantsMin, respRates, respAmps, pathBase, ...
                        x, yLimitsPleth, yLimitsRespRate, yLimitsRespAmp), ...
                    xLimitsMinAll);
else
    handles = struct;
end

%% Output results
respParams.traceFileName = traceFileName;
respParams.pathBase = pathBase;
respParams.plethWindowSeconds = plethWindowSeconds;
respParams.resolutionSeconds = resolutionSeconds;

respData.endPointsWindows = endPointsWindows;
respData.timeInstantsSec = timeInstantsSec;
respData.respRates = respRates;
respData.respAmps = respAmps;

% Place traces in a table
respTable = table(timeInstantsSec, respRates, respAmps);

% Write to a csv file
respTablePath = sprintf('%s_resp_table.csv', pathBase);
writetable(respTable, respTablePath);

% Save outputs to a .mat file
respMatPath = sprintf('%s_resp_data.mat', pathBase);
save(respMatPath, 'respParams', 'respData', '-v7.3');

% Output variably
varargout{1} = respParams;
if nargout >= 2
    varargout{2} = respData;
end
if nargout >= 3
    varargout{3} = handles;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [respRateRpm, respAmp] = ...
                parse_pleth_trace_helper (vector, endPoints, ...
                                            siSeconds, pathBase, ...
                                            plotFlag, respRateRangeRpm)
%% Parses one pleth trace window

%% Hard-coded parameters
smoothWindowRpm = 1;
bandWidthRpm = 2;
xLimitsRpm = [0, 150];

% Convert to Hz
[smoothWindowHz, bandWidthHz, respRateRangeHz] = ...
    argfun(@(x) convert_units(x, 'rpm', 'Hz'), ...
            smoothWindowRpm, bandWidthRpm, respRateRangeRpm);

% Extract the file base
fileBase = extract_fileparts(pathBase, 'base');

% Create a title base
titleBase = replace(fileBase, '_', '\_');

% Extract the piece
piece = extract_subvectors(vector, 'EndPoints', endPoints);

% Subtract the piece by its mean
pieceCentered = piece - mean(piece);

% Compute the number of samples
nSamples = numel(pieceCentered);

% Compute the power spectral density of the piece
% TODO: Choose peak frequencies within the range 30-150 /min
[psdParams, psdData] = ...
    parse_psd(pieceCentered, 'SamplingFrequencyHz', 1/siSeconds, ...
                                'FilterWindowHz', smoothWindowHz, ...
                                'PeakRangeHz', respRateRangeHz);

% Extract peak frequency as respiratory rate
respRateHz = psdParams.peakFrequency;

freqVecHz = psdData.freqVec;
psd = psdData.psd;
psdFiltered = psdData.psdFiltered;
psdSmoothed = psdData.psdSmoothed;

% Filter the trace about the peak frequency
if isnan(respRateHz)
    pieceFiltered = pieceCentered;
else
    % cutoffFreq = respRateHz + bandWidthHz/2;
    cutoffFreq = respRateHz + [-1, 1] * bandWidthHz/2;
    pieceFiltered = freqfilter(pieceCentered, cutoffFreq, siSeconds);
end

% endPointsMiddle = [ceil(nSamples/4); floor(nSamples * 3/4)];
% endPointsMiddle = [1; nSamples];
% pieceMiddle = extract_subvectors(pieceFiltered, 'EndPoints', endPointsMiddle);

% Create a time vector
tVec = siSeconds * (1:nSamples)';
% tVecMiddle = extract_subvectors(tVec, 'EndPoints', endPointsMiddle);

% Compute the respiratory amplitude
% respAmp = max(pieceMiddle) - min(pieceMiddle);
% TODO: Use compute_time_average.m
respAmp = trapz(tVec, abs(pieceCentered)) / (tVec(end) - tVec(1));

% Convert to rpm
respRateRpm = convert_units(respRateHz, 'Hz', 'rpm');

% Plot stuff
if plotFlag
    % Create strings
    respRateStr = sprintf('Resp rate = %.0f / min', respRateRpm);
    respAmpStr = sprintf('Resp amp = %.3g Volts', respAmp);

    % Convert to rpm
    freqVecRpm = convert_units(freqVecHz, 'Hz', 'rpm');

    [fig, ax] = create_subplots(2, 1, 'AlwaysNew', true, ...
                                'FigExpansion', [1, 1]);
    subplot(ax(1));
    plot_traces(tVec, pieceCentered, 'Color', 'k', 'XLabel', 'Time (s)', ...
                'YLabel', 'Pleth (V)', 'FigTitle', titleBase);
    plot_traces(tVec, pieceFiltered, 'Color', 'b', 'PlotOnly', true);
    plot_horizontal_line(respAmp, 'Color', 'g');
    % plot_traces(tVecMiddle, pieceMiddle, 'Color', 'b', 'PlotOnly', true);

    subplot(ax(2));
    plot_traces(freqVecRpm, psd, 'Color', 'k', 'XLabel', 'Frequency (rpm)', ...
                'YLabel', 'Power (V^2/Hz)', 'FigTitle', 'Power spectrum');
    plot_traces(freqVecRpm, psdFiltered, 'Color', 'b', 'PlotOnly', true);
    plot_traces(freqVecRpm, psdSmoothed, 'Color', 'r', 'PlotOnly', true);
    text(0.05, 0.95, respRateStr, 'Units', 'normalized');
    text(0.05, 0.85, respAmpStr, 'Units', 'normalized');
    xlim(xLimitsRpm);

    save_all_figtypes(fig, pathBase, {'png'});
    close(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_resp_traces (timeVecsMin, vector, timeInstantsMin, ...
                                respRates, respAmps, pathBase, xLimitsMin, ...
                                yLimitsPleth, yLimitsRespRate, yLimitsRespAmp)

% Create figure path
figSuffix = sprintf('_pleth_traces_%g-%gmin.csv', ...
            xLimitsMin(1), xLimitsMin(2));
figPath = strcat(pathBase, figSuffix);

% Create figure
[fig, ax] = create_subplots(3, 1, 'AlwaysNew', true, ...
                            'FigExpansion', [1, 1]);

% Link x axis
linkaxes(ax, 'x');

% Plot pleth trace
subplot(ax(1))
plot_traces(timeVecsMin, vector, 'Color', 'k', ...
            'XLimits', xLimitsMin, 'YLimits', yLimitsPleth, ...
            'XLabel', 'Time (min)', 'YLabel', 'Pleth (V)', ...
            'FigTitle', 'suppress');

% Plot resp rate
subplot(ax(2))
plot_traces(timeInstantsMin, respRates, 'Color', 'k', ...
            'XLimits', xLimitsMin, 'YLimits', yLimitsRespRate, ...
            'XLabel', 'Time (min)', 'YLabel', 'Resp rate (rpm)', ...
            'FigTitle', 'suppress');

% Plot resp amp
subplot(ax(3))
plot_traces(timeInstantsMin, respAmps, 'Color', 'k', ...
            'XLimits', xLimitsMin, 'YLimits', yLimitsRespAmp, ...
            'XLabel', 'Time (min)', 'YLabel', 'Resp amp (V)', ...
            'FigTitle', 'suppress');

% Extract the file base
fileBase = extract_fileparts(pathBase, 'base');

% Plot title
suptitle(sprintf('%s: %g-%g min', replace(fileBase, '_', '\_'), ...
                xLimitsMin(1), xLimitsMin(2)));

% Save figure
save_all_figtypes(fig, figPath, {'png'});

handles.fig = fig;
handles.ax = ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
