function plot_FI (varargin)
%% From a current injection protocol, detect spikes for each sweep and make an F-I plot
% Usage: plot_FI (vVecs (opt), iVecs (opt), varargin)
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
% parse_abf
% check_dir
%       cd/construct_and_check_abfpath.m
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
% 2017-05-01 - BT - Converted rangeCI to use identify_CI_protocol.m
% 2017-06-16 - AL - Updated 
% 2018-01-24 - Added isdeployed
% TODO: Change the arguments structure
% TODO: Improve documentation and code structure
% TODO: Use parse_abf.m
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
                ['iVecs must be either a numeric array', ...
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
        siUs = parseParams.siUs;
        tVec = parsedData.tVec;
        iVecs = parsedData.iVecs;
        vVecs = parsedData.vVecs;
    end
else
    % TODO: Organize code
    [tVec, siMs] = match_time_info (tVec, siMs)
    siUs = siMs * 1000;
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

%% Detect spikes
% Determine whether each sample point is a spike for each sweep
isSpike = arrayfun(@(x) detect_spikes_one_sweep(vVecs(:, x)), ...
                    transpose(1:nSweeps), 'UniformOutput', false);

% Put together into a matrix
isSpike = force_matrix(isSpike);

% TODO: Organize code below
%% Compute spike frequencies
% Compute spike frequency for each sweep
spikeFreq = zeros(1,nSweeps);
[~, rangeCI] = identify_CI_protocol(iVecs, siUs);
parfor iSwp = 1:nSweeps
    cdata = vVecs(:, iSwp);
    indSpikes = find(isSpike(:, iSwp));
    indSpikesWithinRange = indSpikes(indSpikes > rangeCI(1) & ...
                                        indSpikes < rangeCI(2));
    if length(indSpikesWithinRange) > 1
        spikeFreq(iSwp) = (length(indSpikesWithinRange) - 1) / ...
                            ((tVec(rangeCI(2)) - tVec(rangeCI(1)))) * 1000;
    else
        spikeFreq(iSwp) = 0;
    end
end

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
    line([tVec(rangeCI(1)) tVec(rangeCI(1))], [y(1) y(2)]);
    line([tVec(rangeCI(2)) tVec(rangeCI(2))], [y(1) y(2)]);
    plot(tVec(indSpikes(indSpikes > rangeCI(1) & indSpikes < rangeCI(2))), ...
          vVecs(indSpikes(indSpikes > rangeCI(1) & indSpikes < rangeCI(2))), 'xg');
    if length(indSpikes) > 1
        text(0.6, 0.85, ['Spike Frequency: ' num2str(spikeFreq(iSwp)) ' Hz'], 'Units', 'normalized');
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
currents = zeros(1, nSweeps);
parfor iSwp = 1:nSweeps
    currents(iSwp) = mean(iVecs(rangeCI(1):rangeCI(2), iSwp));
end
h = figure(999);
set(h, 'Visible', 'off');
clf(h);
plot(currents, spikeFreq);
title(['F-I plot for ', strrep(fileBases, '_', '\_')]);
xlabel('Current Injected (pA)');
ylabel('Spike Frequency (Hz)');
saveas(h, fullfile(outFolder, figNameFI), 'png');
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isSpike = detect_spikes_one_sweep(vVec)
%% Determines whether each sample point is a spike
%   Note: Criteria for a spike: 
%   (1) Must be a local maximum 10 mV higher than the previous local minimum
%   (2) Must be 5 mV higher than the minimum value between the spike 
%           and the following spike

%% Hard-coded parameters
minAmpBefore = 10;     % in mV
minAmpAfter = 5;       % in mV

%% Preparation
% Count the number of samples
nSamples = numel(vVec);

% Initialize array
isSpike = false(nSamples, 1);

%% Do the job
% No spikes if nSamples is less than 3
if nSamples <= 2
    return
end

% Determine whether each data point is a "local maximum" or a "local minimum" 
%   in the general sense
%   TODO: Use findpeaks instead?
isLocalMaximum = [false; diff(vVec(1:end-1)) > 0 & diff(vVec(2:end)) <= 0; ...
                    false];
isLocalMinimum = [false; diff(vVec(1:end-1)) < 0 & diff(vVec(2:end)) >= 0; ...
                    false];

% Finds all peaks that satisfy criterion 1
for idx = 2:nSamples-1
    % Only check local maxima
    if isLocalMaximum(idx)
        % Find the index of the previous local minimum
        idxPrevLocalMinimum = idx - 1;
        while ~isLocalMinimum(idxPrevLocalMinimum) && idxPrevLocalMinimum > 1
            idxPrevLocalMinimum = idxPrevLocalMinimum - 1;
        end

        % Make sure a previous local minimum was found
        if idxPrevLocalMinimum > 1
            % Determine whether criterion 1 is satisfied
            isSpike(idx) = vVec(idx) - vVec(idxPrevLocalMinimum) > minAmpBefore;
        end
    end
end

% Finds all peaks that satisfy criterion 2 by filtering out those 
%   that pass criterion 1 in reverse order
for idx = nSamples-1:-1:2
    % Only check peaks that satisfy criterion 1
    if isSpike(idx)
        % Find the index of the following spike
        idxNextSpike = idx + 1;
        while ~isSpike(idxNextSpike) && idxNextSpike < nSamples
            idxNextSpike = idxNextSpike + 1;
        end

        % Determine whether criterion 2 is satisfied
        isSpike(idx) = vVec(idx) - min(vVec(idx:idxNextSpike)) > minAmpAfter;
    end
end

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