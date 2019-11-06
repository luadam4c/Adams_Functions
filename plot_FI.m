function plot_FI (varargin)
%% From a current injection protocol, detect spikes for each sweep and make an F-I plot
% Usage: plot_FI (vVecs (opt), iVecs (opt), tVec (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       TODO
%
% Arguments:
%       vVecs       - (opt) voltage vectors
%                   must be a numeric array or a cell array of numeric vectors
%                   default == read from .abf file(s)
%       iVecs       - (opt) current vectors
%                   must be a numeric array or a cell array of numeric vectors
%                   default == read from .abf file(s)
%       tVec        - (opt) time vector
%                   must be a numeric array or a cell array of numeric vectors
%                   default == read from .abf file(s)
%       varargin    - 'Directory': the name of the directory containing 
%                                   the abf files, e.g. '20161216'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileNames': names of .abf files to detect
%                   must be a cell array of character arrays or strings
%                   default == detect from pwd
%                   - 'OutFolder': the name of the directory that 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == same as directory
%
% Requires:
% iscellnumeric
% all_files
% extract_fileparts
% extract_subvectors
% match_time_info
% parse_abf
% check_dir
%       cd/construct_and_check_abfpath.m
%       cd/detect_spikes_current_clamp.m
%       cd/identify_channels.m
%       cd/identify_CI_protocol.m
%       /home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%
% Used by:
%       cd/plot_all_abfs.m

% File History:
% 2017-04-04 - Adapted from analyzeCI_20140730
% 2017-04-04 - Used dirr.m to find all abf files
% 2017-04-04 - Updated file names
% 2017-04-06 - BT - Computed and plotted spike frequency
% 2017-04-11 - Moved a copy of identify_channels.m to /home/Matlab/Adams_Functions/
% 2017-04-11 - Now uses fileName directly and calls each file from plot_all_abfs_dir.m
% 2017-04-11 - Cleaned code
% 2017-04-11 - Changed the color map to lines
% 2017-04-13 - BT - Marked on plot spike frequency time interval
% 2017-04-13 - Added alldata, siUs as optional arguments
% 2017-05-01 - BT - Converted endPointsPulse to use identify_CI_protocol.m
% 2017-06-16 - AL - Updated 
% 2018-01-24 - Added isdeployed
% 2019-11-05 - Changed the arguments structure
% 2019-11-05 - Improved documentation and code structure
% 2019-11-05 - Now uses parse_abf.m
% TODO: The injection period isn't right in some cases?
%       See /media/ashleyX/Recordings/20180910/2018_09_10_A_0000_traces

%% Hard-coded parameters
verbose = true;

%% Default values for optional arguments
vVecsDefault = [];              % set later
iVecsDefault = [];              % set later
tVecDefault = [];               % set later
directoryDefault = pwd;         % look for .abf files in 
                                %   the present working directory by default
fileNamesDefault = {};          % detect from pwd by default
outFolderDefault = '';          % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'vVecs', vVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addOptional(iP, 'iVecs', iVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['iVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addOptional(iP, 'tVec', tVecDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
vVecs = iP.Results.vVecs;
iVecs = iP.Results.iVecs;
tVec = iP.Results.tVec;
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
outFolder = iP.Results.OutFolder;

%% Preparation
% Prepare
if isempty(vVecs) || isempty(iVecs) || isempty(tVec)
    % Decide on the file(s) to use
    if isempty(fileNames)
        % Find all .abf files in the directory
        [~, fileNames] = all_files('Directory', directory, ...
                                    'Extension', '.abf', ...
                                    'Verbose', verbose);
    else
        fileNames = construct_and_check_abfpath(fileNames, ...
                                                'Directory', directory);
    end

    % If no files found, return
    if isempty(fileNames)
        fprintf('Type ''help %s'' for usage\n', mfilename);
        return
    end

    % Extract the voltage and current data
    % TODO: What if data is from multiple files
    if iscell(fileNames)
        disp('Not implemented yet!');
        return
    else
        % Parse the .abf file
        [parseParams, parsedData] = parse_abf(fileNames);

        % Extract the time, current and voltage vectors
        siMs = parseParams.siMs;
        tVec = parsedData.tVec;
        iVecs = parsedData.iVecs;
        vVecs = parsedData.vVecs;
    end
else
    % Compute the sampling interval in ms
    [tVec, siMs] = match_time_info(tVec);
end

if ~isempty(fileNames)
    % Extract the directories and file bases
    fileDirs = extract_fileparts(fileNames, 'directory');
    fileBases = extract_fileparts(fileNames, 'base');
else
    fileDirs = directory;
    fileBases = 'some_file';
end

% Set default output folder(s)
if isempty(outFolder)
    outFolder = fullfile(fileDirs, strcat(fileBases, '_traces'));
end

% Create outFolder if not already exists
check_dir(outFolder);

% Count the number of sweeps
nSweeps = size(vVecs, 2);

% Count the number of samples
nSamples = size(vVecs, 1);

%% Parse the current injection protocol
% Identify the (common) current pulse endpoints
[~, endPointsPulse] = identify_CI_protocol(iVecs, siMs);

% Extract the current traces during the injection
iVecsInjected = extract_subvectors(iVecs, 'EndPoints', endPointsPulse);

% Compute current injection values 
currentValues = transpose(mean(iVecsInjected, 1));

%% Detect spikes
[spikesParams, spikesData] = detect_spikes_current_clamp(vVecs);
isSpike = spikesData.isSpike;
idxSpikes = spikesData.idxSpikes;

%% Compute spike frequencies during the injection
% Force as a column cell array
idxSpikes = force_column_cell(idxSpikes);

% Compute the spikes within range
idxSpikesWithinRange = ...
    cellfun(@(x) x(x > endPointsPulse(1) & x < endPointsPulse(2)), ...
            idxSpikes, 'UniformOutput', false);

% Compute spike frequency for each sweep
freqValues = compute_spike_frequency(idxSpikesWithinRange, siMs);

%% Plots
% Plot all sweeps together
h = figure(nSweeps + 1);
set(h, 'Visible', 'off');
clf(h);
figName = [fileBases, '_all'];
cm = colormap(lines);
for iSwp = 1:nSweeps
    plot(tVec, vVecs(:, iSwp), 'Color', cm(mod(iSwp, size(cm, 1)) + 1, :), ...
        'Displayname', ['Sweep #', num2str(iSwp)]);
    hold on;
end
title(strrep(figName, '_', '\_'));
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
legend('Location', 'northeast');
saveas(h, fullfile(outFolder, figName), 'png');
close(h);

% Plot each sweep individually with detected spikes and computed spike frequency    
parfor iSwp = 1:nSweeps
    indSpikes = find(isSpike(:, iSwp));
    figName = [fileBases, '_sweep', num2str(iSwp)];

    h = figure(iSwp);
    set(h, 'Visible', 'off');
    clf(h);
    plot(tVec, vVecs(:, iSwp), 'k')
    hold on;
    plot(tVec(indSpikes), vVecs(indSpikes, iSwp), 'xr');
    y = ylim;
    line([tVec(endPointsPulse(1)) tVec(endPointsPulse(1))], [y(1) y(2)]);
    line([tVec(endPointsPulse(2)) tVec(endPointsPulse(2))], [y(1) y(2)]);
    plot(tVec(indSpikes(indSpikes > endPointsPulse(1) & indSpikes < endPointsPulse(2))), ...
          vVecs(indSpikes(indSpikes > endPointsPulse(1) & indSpikes < endPointsPulse(2))), 'xg');
    if length(indSpikes) > 1
        text(0.6, 0.85, ['Spike Frequency: ' num2str(freqValues(iSwp)) ' Hz'], 'Units', 'normalized');
    else
        text(0.6, 0.85, ['Spike Frequency: 0 Hz'], 'Units', 'normalized');
    end
    title(strrep(figName, '_', '\_'));
    xlabel('Time (ms)')
    ylabel('Membrane Potential (mV)')
    saveas(h, fullfile(outFolder, figName), 'png');
    close(h);
end

% Plot spike frequency over current injected (F-I plot)
figNameFI = [fileBases, '_FI'];
h = figure(999);
set(h, 'Visible', 'off');
clf(h);
plot(currentValues, freqValues);
title(['F-I plot for ', strrep(fileBases, '_', '\_')]);
xlabel('Current Injected (pA)');
ylabel('Spike Frequency (Hz)');
saveas(h, fullfile(outFolder, figNameFI), 'png');
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spikeFreqs = compute_spike_frequency(indSpikes, siMs)
%% Computes the spike frequency for sets of spike indices

% Compute the spike frequency for each set of spike indices
if iscell(indSpikes)
    spikeFreqs = cellfun(@(x) compute_spike_frequency_helper(x, siMs), ...
                        indSpikes);
else
    spikeFreqs = compute_spike_frequency_helper(indSpikes, siMs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spikeFreqHz = compute_spike_frequency_helper(indSpikes, siMs)
%% Computes the spike frequency for one set of spike indices

% Count the number of spikes
nSpikes = numel(indSpikes);

% If less than two spikes, spike frequency is zero
if nSpikes < 2
    spikeFreq = 0;
    return
end

% Otherwise, extract the first and last index
idxFirst = indSpikes(1);
idxLast = indSpikes(end);

% Compute the number of samples between the first and last spike
nSamplesBetweenFirstAndLast = idxLast - idxFirst;

% Compute the time in seconds between the first and last spike
timeDiffSeconds = nSamplesBetweenFirstAndLast * siMs / MS_PER_S;

% Compute the average spike frequency in Hz
spikeFreqHz = (nSpikes - 1) / timeDiffSeconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

isLocalMaximum = false(nSamples, nSweeps);
isLocalMinimum = false(nSamples, nSweeps);
for iSample = 2:nSamples-1
    % Find local maxima
    isLocalMaximum(iSample, :) = ...
        alldata(iSample, :) > alldata(iSample-1, :) & ...
        alldata(iSample, :) >= alldata(iSample+1, :);

    % Find local minima
    isLocalMinimum(iSample, :) = ...
        alldata(iSample, :) < alldata(iSample-1, :) & ...
        alldata(iSample, :) <= alldata(iSample+1, :);
end

isLocalMaximum = ...
    [false(1, nSweeps); ...
    diff(vVecs(1:end-1, :)) > 0 && diff(vVecs(2:end, :)) <= 0; ...
    false(1, nSweeps)];
isLocalMinimum = ...
    [false(1, nSweeps); ...
    diff(vVecs(1:end-1, :)) < 0 && diff(vVecs(2:end, :)) >= 0; ...
    false(1, nSweeps)];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%