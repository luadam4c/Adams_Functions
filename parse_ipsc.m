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
%                   - 'FileBase': base of filename (without extension)
%                   must be a string scalar or a character vector
%                   default == 'unnamed'
%                   - 'StartWindowMs': window (ms) in which IPSC start would lie
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
%       cd/check_subdir.m
%       cd/compute_combined_trace.m
%       cd/create_labels_from_numbers.m
%       cd/extract_common_prefix.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/find_window_endpoints.m
%       cd/match_time_info.m
%       cd/movingaveragefilter.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-11-13 Adapted from find_IPSC_peak.m and find_istart.m
% 2019-11-14 Now uses find_closest.m
% 

%% Hard-coded parameters
% Subdirectories in outFolder for placing figures
outSubDirs = {'ipsc_start', 'ipsc_peak'};

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
startWindowMsDefault = [];      % set later
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
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'StartWindowMs', startWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PeakWindowMs', peakWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
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
startWindowMs = iP.Results.StartWindowMs;
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
if isempty(fileBase)
    fileBase = create_labels_from_numbers(transpose(1:nVectors), ...
                                            'Prefix', 'unnamed_');
end

% Extract common prefix
commonPrefix = extract_common_prefix(fileBase);

% Compute sampling interval(s) and create time vector(s)
if isempty(siMs) && isempty(tVecs)
    error('One of siMs and tVecs must be provided!');
else
    [tVecs, siMs, nSamples] = ...
        match_time_info(tVecs, siMs, nSamples, 'TimeUnits', 'ms')
end

% Compute the base and end of the time vector(s)
tBase = extract_elements(tVecs, 'first') - siMs;
tEnd = extract_elements(tVecs, 'last');

% For verbose
if ~isempty(commonPrefix)
    messageSuffix = ['for ', commonPrefix, ' '];
else
    messageSuffix = '';
end

% Display message
if verbose
    fprintf('ANALYZING IPSC traces for %s ...\n', commonPrefix);
    fprintf('Number of sweeps == %d\n', nVectors);
end

%% Deal with IPSC start
% Decide on the window to look for IPSC start
if isempty(startWindowMs)
    startWindowMs = [mean(tBase); mean(tEnd)]
end

% Detect stimulation start time if not provided
if isempty(stimStartMs)
    % Display message
    if verbose
        fprintf('FINDING common time of IPSC start for %s ...\n', ...
                commonPrefix);
    end

    % Compute the standard deviation trace over all vectors
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
    tVec = tVecs{1};

    % Find end points corresponding to startWindowMs
    endPointsToFindStart = find_window_endpoints(startWindowMs, tVec);

    % Restrict to startWindowMs
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
    peakWindowMs = transpose([stimStartMs, tEnd]);
end

% Find the end points for peakWindowMs
endPointsToFindPeak = find_window_endpoints(peakWindowMs, tVecs);

% Extract the regions of interest
iVecsPart = extract_subvectors(iVecs, 'EndPoints', endPointsToFindPeak);

% Find the minimum in each region
[idxPeakRel, peakAmplitude] = cellfun(@min, iVecsPart);

% Compute the beginning index of the end points for peakWindowMs
idxBegins = extract_elements(endPointsToFindPeak, 'first');

% Compute the original index of the IPSC peak
idxPeak = (idxBegins - 1) + idxPeakRel;

% Compute IPSC peak time in ms
peakTimeMs = extract_elements(tVecs, 'specific', 'Index', idxPeak);

% Compute IPSC peak delay in ms
peakDelayMs = peakTimeMs - stimStartMs;

%% Output results in tables
% TODO

%% Plot current trace
if plotflag
    % Check if needed outSubDirs exist in outfolder
    check_subdir(outfolder, outSubDirs);

    % Plot current traces
    fig = set_figure_properties('Name', 'Current analysis', ...
                            'AlwaysNew', true);
    if ~isempty(vVec0s)
        subplot(2,1,1);
    end
    for iSwp = 1:nVectors                   % Plot each current trace
        plot(tVecs(ind), iVecs(ind, iSwp)); hold on; 
    end
    plot(tVecs(ind), iStdVec(ind), 'LineWidth', 2, 'Color', 'r');
    plot(tVecs(ind), iStdVecSmooth(ind), 'LineWidth', 2, 'Color', 'g');
    plot(tVecs(ind), diffIStdVecsSmooth(ind), 'LineWidth', 2, 'Color', 'b');
    if ~isempty(idxRelFirstChange)
        for iSwp = 1:nVectors
            plot(tVecs(idxStimStart), iVecs(idxStimStart, iSwp), ...
                'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
        end
    end
    title(['Current analysis for ', strrep(filebase, '_', '\_')]);
    xlabel('Time (ms)')
    ylabel('Current (pA)')
    xlim(startWindowMs);
%   ylim([-30 30]);             
    if  ~isempty(vVec0s)
        subplot(2,1,2);
        for iSwp = 1:nVectors       % Plot each voltage trace
            plot(tVecs(ind), vVec0s(ind, iSwp)); hold on; 
        end
        if ~isempty(idxRelFirstChange)
            for iSwp = 1:nVectors
                plot(tVecs(idxStimStart), vVec0s(idxStimStart, iSwp), ...
                    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
            end
        end
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
        xlim(startWindowMs);
%       ylim([-90 40]);             
    end
    figname = fullfile(outfolder, directories{1}, [filebase, '_istart', '.png']);
    saveas(fig, figname);
    close(fig);

    % Plot current trace
    fig = figure(4000);
    set(fig, 'Visible', 'off');
    set(fig, 'Name', 'IPSC peak amplitude analysis');
    clf(fig);
    for iSwp = 1:nVectors            % Plot each current trace and mark peak amplitude
        plot(tVec0(ind), iVecs(ind, iSwp)); hold on; 
        plot(tVec0(idxPeak(iSwp)), peakAmplitude(iSwp), 'Marker', 'x', 'MarkerSize', 12);
    end
    xlabel('Time (ms)')
    ylabel('Current (pA)')
    title(['IPSC peak amplitude analysis for ', strrep(filebase, '_', '\_')]);
    figname = fullfile(outfolder, outSubDirs{2}, [filebase, '_IPSCpeak', '.png']);
    saveas(fig, figname);
    close(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

IPSC_ind = round(ipscpwin(1)/siMs);     % Assume no IPSC offset % NOT ROBUST
ind = IPSC_ind:round(ipscpwin(2)/siMs);     % indices of interest   % NOT ROBUST

    [peakAmplitude, itemp1] = min(iVec0s_part1);    % peakAmplitude is the amplitude of ipeak
    idxPeak = (ind(1) - 1) + itemp1;      % index of ipeak relative to iVecs
    peakTimeMs = tVec0(idxPeak);          % time of ipeak in ms

    close(fig);

if isempty(startWindowMs)
    if ~isempty(stimStartMs)
        stimStartMsRow = force_row_vector(stimStartMs);
        startWindowMs = ones(2, 1) * stimStartMsRow;
    else
        startWindowMs = transpose([tBase, tEnd]);
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
