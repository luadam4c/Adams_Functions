function varargout = parse_pleth_trace (vectors, siSeconds, varargin)
%% Parses pleth traces
% Usage: [parsedParams, parsedData] = parse_pleth_trace (vectors, siSeconds, varargin)
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
%       timeVecs = create_time_vectors('TimeUnits', 's', 'SamplingIntervalSeconds', siSeconds, 'Vectors', plethVec);
%       [parsedParams, parsedData] = parse_pleth_trace(plethVec, siSeconds, 'TraceFileName', spike2MatPath);
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       vectors     - pleth trace(s)
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       siSeconds   - sampling interval in seconds
%                   must be a positive vector
%       varargin    - 'TraceFileName': Name of the corresponding trace file(s)
%                   must be empty, a character vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - 'FileBase': file base for output files
%                   must be a string scalar or a character vector
%                   default == match the trace file name
%                   - 'PlethWindowSeconds' - pleth window for computing
%                                               breathing rate in seconds
%                   must be a numeric scalar
%                   default = 60 seconds
%                   - 'ResolutionSeconds' - resolution in seconds
%                   must be a numeric scalar
%                   default = 60 seconds
%                   - Any other parameter-value pair for TODO
%
% Requires:
%       cd/array_fun.m
%       cd/compute_running_windows.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_vectors.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/freqfilter.m
%       cd/isemptycell.m
%       cd/match_row_count.m
%       cd/match_format_vector_sets.m
%       cd/parse_psd.m
%
% Used by:
%       cd/parse_spike2_mat.m

% File History:
% 2020-08-12 Modified from parse_laser_trace.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
traceFileNameDefault = '';      % set later
plethWindowSecondsDefault = 60;
resolutionSecondsDefault = 60;
fileBaseDefault = '';

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
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siSeconds', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TraceFileName', traceFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlethWindowSeconds', plethWindowSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ResolutionSeconds', resolutionSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, vectors, siSeconds, varargin{:});
traceFileName = iP.Results.TraceFileName;
fileBase = iP.Results.FileBase;
plethWindowSeconds = iP.Results.PlethWindowSeconds;
resolutionSeconds = iP.Results.ResolutionSeconds;

% Keep unmatched arguments for the TODO function
otherArguments = iP.Unmatched;

%% Preparation
% Force vectors to be column vectors
vectors = force_column_vector(vectors);

% Count the number of samples
nSamples = count_samples(vectors);

% Count the number of vectors
nVectors = count_vectors(vectors);

% Create time vectors in seconds
timeVecs = create_time_vectors(nSamples, 'TimeUnits', 's', ...
                    'SamplingIntervalSeconds', siSeconds, 'Vectors', vectors);

% Match the row count
siSeconds = match_row_count(siSeconds, nVectors);

% Match the number of file names to the number of vectors
[traceFileName, vectors] = match_format_vector_sets(traceFileName, vectors);

% Make sure fileBase is in agreement with nVectors
if isempty(fileBase)
    if isempty(traceFileName)
        fileBase = create_time_stamp;
    elseif isemptycell(traceFileName)
        fileBase = create_labels_from_numbers(1:nVectors, ...
                                'Prefix', strcat(create_time_stamp, '_'));
    else
        % Extract from the trace file name
        fileBase = extract_fileparts(traceFileName, 'pathbase');
    end
else
    if ischar(fileBase) && nVectors > 1 || ...
            iscell(fileBase) && numel(fileBase) ~= nVectors
        fprintf('Number of file bases must match the number of vectors!\n');
        varargout{1} = table.empty;
        varargout{2} = table.empty;
        return
    end
end

%% Do the job
% Compute running windows
[endPointsWindows, timeInstants] = ...
    compute_running_windows(timeVecs, plethWindowSeconds, 'IsRegular', true, ...
                                'Resolution', resolutionSeconds);

% Count the number of windows
nWindows = size(endPointsWindows, 2);

% Create new file bases
fileBases = create_labels_from_numbers(1:nWindows, ...
                            'Prefix', strcat(fileBase, '_pleth_window'));

% Create window numbers
windowNumbers = (1:nWindows)';

% Decide whether to plot each window
plotFlags = mod(windowNumbers, 1000) == 0;

% Compute respiratory rates and amplitudes for each running window
[respRates, respAmps] = ...
    array_fun(@(x, y) parse_pleth_trace_helper(vectors, endPointsWindows(:, x), ...
                            siSeconds, fileBases{x}, y), ...
        	windowNumbers, plotFlags);

%% Plot results
% Decide on x limits to plot
xLimits = [0, 2400];

% Create figure path
figSuffix = sprintf('_pleth_traces_%g-%gsec.csv', xLimits(1), xLimits(2));
figPath = strcat(fileBase, figSuffix);

% Create figure
[fig, ax] = create_subplots(3, 1, 'AlwaysNew', true, 'FigExpansion', [1, 1]);

% Link x axis
linkaxes(ax, 'x');

% Plot pleth trace
subplot(ax(1))
plot_traces(timeVecs, vectors, 'Color', 'k', ...
            'XLabel', 'Time (s)', 'YLabel', 'Pleth (V)', ...
            'FigTitle', 'suppress');

% Plot resp rate
subplot(ax(2))
plot_traces(timeInstants, respRates, 'Color', 'k', ...
            'XLabel', 'Time (s)', 'YLabel', 'Resp rate (Hz)', ...
            'FigTitle', 'suppress');

% Plot resp amp
subplot(ax(3))
plot_traces(timeInstants, respAmps, 'Color', 'k', ...
            'XLabel', 'Time (s)', 'YLabel', 'Resp amp (V)', ...
            'FigTitle', 'suppress');

% Set x limits
xlim(xLimits);

% Save figure
save_all_figtypes(fig, figPath, {'png'});

%% Output results
parsedParams.plethWindowSeconds = plethWindowSeconds;
parsedParams.resolutionSeconds = resolutionSeconds;

parsedData.endPointsWindows = endPointsWindows;
parsedData.timeInstants = timeInstants;
parsedData.respRates = respRates;
parsedData.respAmps = respAmps;

% Place traces in a table
respTable = table(timeInstants, respRates, respAmps);

% Write to a file
respPath = strcat(fileBase, '_pleth_traces.csv');
writetable(respTable, respPath);

% Output variably
varargout{1} = parsedParams;
if nargout > 1
    varargout{2} = parsedData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [respRate, respAmp] = parse_pleth_trace_helper (vector, endPoints, ...
                                                siSeconds, pathBase, plotFlag)
%% Parses one pleth trace window

%% Hard-coded parameters
smoothWindowHz = 0.1;
bandWidthHz = 1;

% Extract the file base
fileBase = extract_fileparts(pathBase, 'base');

% Create a title base
titleBase = replace(fileBase, '_', '\_');

% Extract the piece
piece = extract_subvectors(vector, 'EndPoints', endPoints);

% Subtract the piece with its mean
pieceCentered = piece - mean(piece);

% Compute the number of samples
nSamples = numel(pieceCentered);

% Compute the power spectral density of the piece
[psdParams, psdData] = ...
    parse_psd(pieceCentered, 'SamplingFrequencyHz', 1/siSeconds, ...
                                'FilterWindowHz', smoothWindowHz);

% Extract peak frequency as respiratory rate
respRate = psdParams.peakFrequency;

freqVec = psdData.freqVec;
psd = psdData.psd;
psdFiltered = psdData.psdFiltered;
psdSmoothed = psdData.psdSmoothed;

% Filter the trace about the peak frequency
% cutoffFreq = respRate + [-1, 1] * bandWidthHz/2;
if isnan(respRate)
    pieceFiltered = pieceCentered;
else
    cutoffFreq = respRate + bandWidthHz/2;
    pieceFiltered = freqfilter(pieceCentered, cutoffFreq, siSeconds);
end

% endPointsMiddle = [ceil(nSamples/4); floor(nSamples * 3/4)];
endPointsMiddle = [1; nSamples];
pieceMiddle = extract_subvectors(pieceFiltered, 'EndPoints', endPointsMiddle);

% Compute the respiratory amplitude
respAmp = max(pieceMiddle) - min(pieceMiddle);

% Plot stuff
if plotFlag
    % Create a time vector
    tVec = siSeconds .* (1:nSamples)';
    tVecMiddle = extract_subvectors(tVec, 'EndPoints', endPointsMiddle);

    [fig, ax] = create_subplots(2, 1, 'AlwaysNew', true, 'FigExpansion', [1, 1]);
    subplot(ax(1));
    plot_traces(tVec, pieceCentered, 'Color', 'k', 'XLabel', 'Time (s)', ...
                'YLabel', 'Pleth (V)', 'FigTitle', titleBase);
    plot_traces(tVecMiddle, pieceMiddle, 'Color', 'b', 'PlotOnly', true);

    subplot(ax(2));
    plot_traces(freqVec, psd, 'Color', 'k', 'XLabel', 'Frequency (Hz)', ...
                'YLabel', 'Power (V^2/Hz)', 'FigTitle', 'Power spectrum');
    plot_traces(freqVec, psdFiltered, 'Color', 'b', 'PlotOnly', true);
    plot_traces(freqVec, psdSmoothed, 'Color', 'r', 'PlotOnly', true);
    xlim([0, 20]);

    save_all_figtypes(fig, pathBase, {'png'});
    close(fig);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Subtract the vector with its mean
vectorCentered = vector - mean(vector);
vectorFiltered = freqfilter(vectorCentered, cutoffFreqs, siSeconds);
pieceFiltered = extract_subvectors(vectorFiltered, 'EndPoints', endPoints);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
