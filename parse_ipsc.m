function varargout = parse_ipsc (iVecs, varargin)
%% Finds time of current start and peak from a an inhibitory post-synaptic current trace (must be negative current)
% Usage: [parsedParams, parsedData] = parse_ipsc (iVecs, siMs (opt), varargin)
% Explanation:
%       TODO
%
% Examples:
%       sweepName = 'C101210_0006_3';
%       matFilesDir = '/media/adamX/m3ha/data_dclamp/take4/matfiles';
%       [data, sweepInfo] = m3ha_import_raw_traces(sweepName, 'Directory', matFilesDir);
%       iVecs = extract_columns(data, 3);
%       siMs = sweepInfo.siMs;
%       [parsedParams, parsedData] = parse_ipsc(iVecs, siMs);
%
% Outputs: 
%       parsedParams    - a table containing the parsed parameters, 
%                           each row corresponding to a vector, with fields:
%                           idxStimStart    - index of IPSC start
%                           stimStartMs     - time of IPSC start in ms
%                           idxPeak         - index of IPSC peak
%                           peakDelayMs     - time delay of IPSC peak in ms
%                           peakTimeMs      - time of IPSC peak in ms
%                           peakAmplitude   - amplitude of IPSC peak in pA
%                       specified as a table
%
% Arguments:
%       iVecs      - original current vector(s) in pA
%                   must be a numeric array
%                       or a cell array of numeric vectors
%       siMs        - (opt) sampling interval in ms
%                           If not provided, 'tVecs' must be provided
%                   must be a positive vector
%                   default == computed from tVecs
%       varargin    - 'Verbose': whether to output parsed results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotFlag': whether to plot traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileBase': base of file name (without extension) 
%                                   corresponding to each vector
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == 'unnamed_1', 'unnamed_2', ...
%                   - 'StimStartWindowMs': window (ms) in which 
%                                               IPSC start would lie
%                   must be a positive scalar
%                   default == [mean(tBase); mean(tEnd)]
%                           Note: mh3a uses [1000, 1100]
%                   - 'StimStartMs': IPSC start time (ms), 
%                           Note: this is the time relative to which 
%                                       the peak delay is computed
%                   must be a positive scalar
%                   default == detected
%                   - 'PeakWindowMs': window (ms) in which IPSC peak would lie
%                   must be a positive scalar
%                   default == transpose([stimStartMs, tEnd])
%                   - 'tVecs': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   default == created from siMs and iVecs
%                   - 'vVecs': original voltage vector(s)
%                   must be a numeric array with same length as vVec0
%                   default == [] (not provided)
%
% Requires:
%       cd/argfun.m
%       cd/check_subdir.m
%       cd/compute_combined_trace.m
%       cd/decide_on_filebases.m
%       cd/extract_common_prefix.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/find_window_endpoints.m
%       cd/match_row_count.m
%       cd/match_time_info.m
%       cd/movingaveragefilter.m
%       cd/plot_traces.m
%       cd/set_figure_properties.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-11-13 Adapted from find_IPSC_peak.m and find_istart.m
% 2019-11-14 Now uses find_closest.m
% 2019-11-15 Now uses match_row_count.m
% TODO: Test this function, especially with PlotFlag true
% 

%% Hard-coded parameters
% Subdirectories in outFolder for placing figures
% outSubDirs = {'IPSCstart', 'IPSCpeak'};
outSubDirs = {'IPSCpeak'};

% Parameters used for data analysis
smoothWindowMs = 0.3;       % width in ms for the moving average filter
                            %    for finding IPSC start
slopeThreshold = 5;         % slope threshold for finding stimStartMs

%% Default values for optional arguments
siMsDefault = [];               % set later
verboseDefault = true;          % print to standard output by default
plotFlagDefault = false;
outFolderDefault = pwd;
fileBaseDefault = {};           % set later
stimStartWindowMsDefault = [];      % set later
stimStartMsDefault = [];        % set later
peakWindowMsDefault = [];       % set later
tVecsDefault = [];              % set later
vVecsDefault = [];              % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'iVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['iVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'siMs', siMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileBase must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'StimStartWindowMs', stimStartWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PeakWindowMs', peakWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVecs', vVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, iVecs, varargin{:});
siMs = iP.Results.siMs;
verbose = iP.Results.Verbose;
plotFlag = iP.Results.PlotFlag;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
stimStartWindowMs = iP.Results.StimStartWindowMs;
stimStartMs = iP.Results.StimStartMs;
peakWindowMs = iP.Results.PeakWindowMs;
tVecs = iP.Results.tVecs;
vVecs = iP.Results.vVecs;

%% Preparation
% Count the number of vectors
nVectors = count_vectors(iVecs);

% Count the number of samples for each vector
nSamples = count_samples(iVecs);

% Create file bases if not provided
fileBase = decide_on_filebases(fileBase, nVectors);

% Extract common prefix
commonPrefix = extract_common_prefix(fileBase);

% Compute sampling interval(s) and create time vector(s)
if isempty(siMs) && isempty(tVecs)
    error('One of siMs and tVecs must be provided!');
else
    [tVecs, siMs] = ...
        match_time_info(tVecs, siMs, nSamples, 'TimeUnits', 'ms');
end

% Compute the base and end of the time vector(s)
tBase = extract_elements(tVecs, 'first') - siMs;
tEnd = extract_elements(tVecs, 'last');

% Display message
if verbose
    fprintf('ANALYZING IPSC traces for %s ...\n', commonPrefix);
    fprintf('Number of sweeps == %d\n', nVectors);
end

%% Deal with IPSC start
% Decide on the window to look for IPSC start
if isempty(stimStartWindowMs)
    stimStartWindowMs = [mean(tBase); mean(tEnd)];
end

% Detect stimulation start time if not provided
% TODO: Fit an exponential? if only one iVec is provided
if isempty(stimStartMs)
    % Save the detection status
    stimStartDetected = true;
    
    % Display message
    if verbose
        fprintf('FINDING common time of IPSC start for %s ...\n', ...
                commonPrefix);
    end

    % If there are more than one vectors, 
    %   compute the standard deviation trace over all vectors
    iStdVec = compute_combined_trace(iVecs, 'std');

    % Compute an average sampling interval
    siMsAveraged = mean(siMs);
    if verbose
        fprintf('Average sampling interval == %g ms\n', siMsAveraged);
    end

    % Smooth the standard deviation trace
    iStdVecSmooth = movingaveragefilter(iStdVec, smoothWindowMs, siMsAveraged);

    % Compute the slope of of the smoothed standard deviation trace
    diffIStdVecsSmooth = diff(iStdVecSmooth) / siMsAveraged;

    % Use the first time vector
    %   TODO: What if time vectors are different?
    if iscell(tVecs)
        tVec = tVecs{1};
    else
        tVec = tVecs;
    end

    % Find end points corresponding to stimStartWindowMs
    endPointsToFindStart = find_window_endpoints(stimStartWindowMs, tVec);

    % Restrict to stimStartWindowMs
    diffIStdVecsSmoothRestricted = extract_subvectors(diffIStdVecsSmooth, ...
                                            'EndPoints', endPointsToFindStart);

    % Find the first time point that the slope of the standard deviation
    %   reaches threshold
    idxRelFirstChange = find(diffIStdVecsSmoothRestricted > slopeThreshold, ...
                                1, 'first');


    % Compute index of IPSC start
    if ~isempty(idxRelFirstChange)
        idxStimStart = endPointsToFindStart(1) + idxRelFirstChange - 1;
    else
        idxStimStart = NaN;
    end

    % Compute IPSC start time in ms
    if isnan(idxStimStart)
        stimStartMs = NaN;
    else
        stimStartMs = tVec(idxStimStart);
    end

    % Display results
    if verbose
        fprintf('Index of current application == %d\n', idxStimStart);
        fprintf('Time of current application == %g ms\n', stimStartMs);
    end
else
    % Stimulation start not detected
    stimStartDetected = false;

    % Use the indices of tVecs with values closest to stimStartMs
    idxStimStart = find_closest(tVecs, stimStartMs);
end

%% Deal with IPSC peak
if verbose
    fprintf('FINDING time of IPSC peak for %s ...\n', commonPrefix);
    fprintf('Sampling interval == %g ms\n', siMs);
end

% Decide on the window to look for IPSC peak
if isempty(peakWindowMs)
    if ~isnan(stimStartMs)
        peakWindowMs = transpose([stimStartMs, tEnd]);
    else
        peakWindowMs = transpose([tBase, tEnd]);
    end
end

% Find the end points for peakWindowMs
endPointsToFindPeak = find_window_endpoints(peakWindowMs, tVecs);

% Extract the regions of interest
iVecsPart = extract_subvectors(iVecs, 'EndPoints', endPointsToFindPeak);

% Find the minimum in each region
[peakAmplitude, idxPeakRel] = cellfun(@min, iVecsPart);

% Compute the beginning index of the end points for peakWindowMs
idxBegins = extract_elements(endPointsToFindPeak, 'first');

% Compute the original index of the IPSC peak
idxPeak = (idxBegins - 1) + idxPeakRel;

% Compute IPSC peak time in ms
peakTimeMs = extract_elements(tVecs, 'specific', 'Index', idxPeak);

% Compute IPSC peak delay in ms
peakDelayMs = peakTimeMs - stimStartMs;

%% Organize results
% Match row counts
[idxStimStart, stimStartMs, idxPeak, ...
        peakDelayMs, peakTimeMs, peakAmplitude] = ...
    argfun(@(x) match_row_count(x, nVectors), ...
            idxStimStart, stimStartMs, idxPeak, ...
            peakDelayMs, peakTimeMs, peakAmplitude);

% Put parameters in a table
parsedParams = table(idxStimStart, stimStartMs, idxPeak, peakDelayMs, ...
                    peakTimeMs, peakAmplitude);

% Put data in a structure
parsedData.tVecs = tVecs;
parsedData.iVecs = iVecs;
if stimStartDetected
    parsedData.iStdVec = iStdVec;
    parsedData.iStdVecSmooth = iStdVecSmooth;
    parsedData.diffIStdVecsSmooth = diffIStdVecsSmooth;
end

%% Output results
varargout{1} = parsedParams;
varargout{2} = parsedData;

%% Plot current traces
% TODO: Not organized and tested:
if plotFlag
    % Check if needed outSubDirs exist in outFolder
    check_subdir(outFolder, outSubDirs);

    % Plot current peak analysis
    figName = fullfile(outFolder, outSubDirs{1}, ...
                        [commonPrefix, '_IPSCpeak', '.png']);
    figTitle = ['IPSC peak amplitude analysis for ', ...
                    replace(commonPrefix, '_', '\_')];
    figHandle = set_figure_properties('AlwaysNew', true, 'Visible', 'off', ...
                                'Name', 'IPSC peak amplitude analysis');
    plot_traces(tVecs, iVecs, 'PlotMode', 'overlapped', ...
                'XLimits', peakWindowMs, ...
                'XLabel', 'Time (ms)', 'YLabel', 'Current (pA)', ...
                'LegendLocation', 'suppress', ...
                'FigHandle', figHandle, 'FigTitle', figTitle, ...
                'LineWidth', 1);
    hold on
    arrayfun(@(x) plot(tVecs{x}(idxPeak(x)), peakAmplitude(x), ...
                        'rx', 'MarkerSize', 8, 'LineWidth', 1), ...
                        transpose(1:nVectors));
    save_all_figtypes(figHandle, figName);
    close(figHandle);

    % Plot current and voltage traces
    % TODO: fix
%{
    % Find the indices to plot
    indToPlot = endPointsToFindPeak(1):endPointsToFindPeak(2);

    figHandle = set_figure_properties('Name', 'Current analysis', ...
                                'AlwaysNew', true);
    if ~isempty(vVecs)
        subplot(2,1,1);
    end
    for iSwp = 1:nVectors                   % Plot each current trace
        plot(tVecs(indToPlot), iVecs(indToPlot, iSwp)); hold on; 
    end
    plot(tVecs(indToPlot), iStdVec(indToPlot), 'LineWidth', 2, 'Color', 'r');
    plot(tVecs(indToPlot), iStdVecSmooth(indToPlot), 'LineWidth', 2, 'Color', 'g');
    plot(tVecs(indToPlot), diffIStdVecsSmooth(indToPlot), 'LineWidth', 2, 'Color', 'b');
    if ~isempty(idxRelFirstChange)
        for iSwp = 1:nVectors
            plot(tVecs(idxStimStart), iVecs(idxStimStart, iSwp), ...
                'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
        end
    end
    title(['Current analysis for ', strrep(fileBase, '_', '\_')]);
    xlabel('Time (ms)')
    ylabel('Current (pA)')
    xlim(stimStartWindowMs);
%   ylim([-30 30]);             
    if  ~isempty(vVecs)
        subplot(2,1,2);
        for iSwp = 1:nVectors       % Plot each voltage trace
            plot(tVecs(indToPlot), vVecs(indToPlot, iSwp)); hold on; 
        end
        if ~isempty(idxRelFirstChange)
            for iSwp = 1:nVectors
                plot(tVecs(idxStimStart), vVecs(idxStimStart, iSwp), ...
                    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
            end
        end
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
        xlim(stimStartWindowMs);
%       ylim([-90 40]);             
    end
    figName = fullfile(outFolder, outSubDirs{1}, [fileBase, '_istart', '.png']);
    saveas(figHandle, figName);
    close(figHandle);
%}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if isempty(stimStartWindowMs)
    if ~isempty(stimStartMs)
        stimStartMsRow = force_row_vector(stimStartMs);
        stimStartWindowMs = ones(2, 1) * stimStartMsRow;
    else
        stimStartWindowMs = transpose([tBase, tEnd]);
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
