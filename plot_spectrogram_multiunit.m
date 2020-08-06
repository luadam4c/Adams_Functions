function handles = plot_spectrogram_multiunit (parsedData, parsedParams, varargin)
%% Plots spectrograms from parsed multiunit data
% Usage: handles = plot_spectrogram_multiunit (parsedData, parsedParams, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - a structure with fields:
%                       TODO
%                   specified as a scalar structure
%
% Arguments:
%       parsedData  - parsed data
%                   must be a table with fields:
%                       tVecs
%                       vVecs
%                       vVecsFilt
%       parsedParams- parsed parameters
%                   must be a table with fields:
%                       stimStartSec
%                       figPathBase
%       varargin    - 'BinWidthSeconds': bin width in seconds
%                   must be a numeric scalar
%                   default == 0.1
%                   - 'PlotStim': whether to plot stimulation time
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SweepNumbers' - the sweep numbers to plot
%                   must be a numeric vector or 'all'
%                   default == 'all'
%                   - 'XLimitsSec': limits of x axis in seconds
%                   must be empty or a 2-element increasing numeric vector
%                   default == [2, 20] seconds
%                   - 'YLimits': limits of y axis in Hz, 
%                               suppress by setting value to empty
%                   must be empty or a 2-element increasing numeric vector
%                   default == [1, 1500] Hz
%                   - 'CLimits': limits of color axis in V s, 
%                               suppress by setting value to empty
%                   must be empty or a 2-element increasing numeric vector
%                   default == [0, 1] V s
%                   - Any other parameter-value pair for plot_spectrogram()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/create_subplots.m
%       cd/compute_spectrogram.m
%       cd/compute_sampling_interval.m
%       cd/convert_units.m
%       cd/find_window_endpoints.m
%       cd/extract_subvectors.m
%       cd/plot_spectrogram.m
%       cd/plot_vertical_line.m
%
% Used by:
%       cd/m3ha_oscillations_analyze.m

% File History:
% 2019-12-02 Adapted from plot_spectrogram_multiunit.m

%% Hard-coded parameters


%% Default values for optional arguments
binWidthSecondsDefault = 0.1;   % 100 ms bins by default
plotStimDefault = true;         % plot stimulation time by default
sweepNumbersDefault = 'all';    % plot all sweeps by default
xLimitsSecDefault = [2, 20];  	% plots 2-20 sec by default
yLimitsDefault = [1, 1500];  	% plots 1-1500 Hz by default
cLimitsDefault = [0, 1];  	    % plots 0-1 V s by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'parsedData', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addRequired(iP, 'parsedParams', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BinWidthSeconds', binWidthSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'PlotStim', plotStimDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SweepNumbers', sweepNumbersDefault);
addParameter(iP, 'XLimitsSec', xLimitsSecDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'CLimits', cLimitsDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, parsedData, parsedParams, varargin{:});
binWidthSeconds = iP.Results.BinWidthSeconds;
plotStim = iP.Results.PlotStim;
sweepNumbers = iP.Results.SweepNumbers;
xLimitsSec = iP.Results.XLimitsSec;
yLimits = iP.Results.YLimits;
cLimits = iP.Results.CLimits;

% Keep unmatched arguments for the plot_spectrogram() function
otherArguments = iP.Unmatched;

%% Preparation
% Restrict to specific sweeps if requested
if isnumeric(sweepNumbers) && ~isempty(sweepNumbers)
    [parsedData, parsedParams] = ...
        argfun(@(x) x(sweepNumbers, :), parsedData, parsedParams);
end

% Extract parameters
tVecsMs = parsedData.tVec;
vVecs = parsedData.vVec;

% Extract parameters
stimStartSec = parsedParams.stimStartSec;
figPathBase = parsedParams.figPathBase;

% Count the number of sweeps
nSweeps = height(parsedParams);

% Convert time vector to seconds
tVecsSec = convert_units(tVecsMs, 'ms', 's');

% Extract the file bases
fileBases = extract_fileparts(figPathBase, 'base');

% Create figure titles
figSubTitles = replace(fileBases, '_', '\_');

%% Plot the spectrograms
% Create subplots
[~, axHandles] = create_subplots(nSweeps, 1, 'FigExpansion', [1, nSweeps/4]);

% Compute and plot each spectrogram
handles = ...
    cellfun(@(a, b, c, d) plot_spectrogram_multiunit_helper(a, b, c, d, ...
                        xLimitsSec, yLimits, cLimits, ...
                        plotStim, stimStartSec, ...
                        binWidthSeconds, otherArguments), ...
            tVecsSec, vVecs, num2cell(axHandles), figSubTitles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_spectrogram_multiunit_helper(tVecSec, vVec, ...
                            axHandle, figSubTitle, ...
                            xLimitsSec, yLimits, cLimits, ...
                            plotStim, stimStartSec, binWidthSeconds, ...
                            otherArguments)

% Go to this subplot
subplot(axHandle);

% Compute the sampling interval in seconds
siSeconds = compute_sampling_interval(tVecSec);

% Deal with x axis limits
if ~isempty(xLimitsSec)
    % The start time is the left x limit
    startTimeSeconds = xLimitsSec(1);
    
    % Find x limits end points
    endPoints = find_window_endpoints(xLimitsSec, tVecSec);

    % Extract the subvector
    vVec = extract_subvectors(vVec, 'EndPoints', endPoints);
else
    % Thd start time is zero
    startTimeSeconds = 0;
end

% Compute the spectrogram in [mV s]
[spectDataMvs, freqHz, timeInstantsSeconds] = ...
    compute_spectrogram(vVec, siSeconds, 'BinWidthSeconds', binWidthSeconds, ...
                        'StartTimeSeconds', startTimeSeconds);

% Convert the spectrogram to [V s]
spectDataVs = convert_units(spectDataMvs, 'mVs', 'Vs');

% Plot the magnitude of the spectrogram in [V s]
handles = plot_spectrogram(spectDataVs, timeInstantsSeconds, freqHz, ...
                            'YLimits', yLimits, 'CLimits', cLimits, ...
                            otherArguments);

% Restrict to x limits
xlim(xLimitsSec);
                        
% Plot stimulation start
if plotStim
    vertLine = plot_vertical_line(mean(stimStartSec), 'Color', 'g', ...
                                    'LineStyle', '--', 'LineWidth', 0.5);
end

% Plot a subtitle
if ~isempty(figSubTitle)
    title(figSubTitle);
end

if plotStim
    handles.vertLine = vertLine;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%