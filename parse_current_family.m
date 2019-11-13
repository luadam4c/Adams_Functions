function [parsedParams, parsedData] = parse_current_family (varargin)
%% From a family of current injections, detect current pulse times and amplitudes, detect spikes for each voltage response, and compute spike frequencies
% Usage: [parsedParams, parsedData] = parse_current_family (vVecs (opt), iVecs (opt), tVec (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       parse_current_family;
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
%                   - 'OutFileBase': file base for output files
%                   must be a string scalar or a character vector
%                   default == match the file name
%
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/check_dir.m
%       cd/compute_spike_frequency.m
%       cd/construct_and_check_abfpath.m
%       cd/count_vectors.m
%       cd/create_labels_from_numbers.m
%       cd/create_subplots.m
%       cd/detect_spikes_current_clamp.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/force_column_cell.m
%       cd/identify_channels.m
%       cd/identify_CI_protocol.m
%       cd/iscellnumeric.m
%       cd/match_row_count.m
%       cd/match_time_info.m
%       cd/parse_abf.m
%       cd/plot_window_boundaries.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%       cd/set_figure_properties.m
%       cd/save_all_figtypes.m
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
% 2017-05-01 - BT - Converted pulseEndPoints to use identify_CI_protocol.m
% 2017-06-16 - AL - Updated 
% 2018-01-24 - Added isdeployed
% 2019-11-05 - Changed the arguments structure
% 2019-11-05 - Improved documentation and code structure
% 2019-11-05 - Now uses parse_abf.m
% TODO: The injection period isn't right in some cases?
%       See /media/ashleyX/Recordings/20180910/2018_09_10_A_0000_traces

%% Hard-coded parameters
% TODO: Make optional arguments
verbose = true;
plotRawData = false; %true;
plotSpikeDetection = false %true;
plotSeparately = false;
plotFI = false; %true;
parsedParamsSuffix = 'current_family_params';
parsedDataSuffix = 'current_family_data';

%% Default values for optional arguments
vVecsDefault = [];              % set later
iVecsDefault = [];              % set later
tVecDefault = [];               % set later
directoryDefault = pwd;         % look for .abf files in 
                                %   the present working directory by default
fileNamesDefault = {};          % detect from pwd by default
outFolderDefault = '';          % set later
outFileBaseDefault = '';        % set later

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
addParameter(iP, 'OutFileBase', outFileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
vVecs = iP.Results.vVecs;
iVecs = iP.Results.iVecs;
tVec = iP.Results.tVec;
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
outFolder = iP.Results.OutFolder;
outFileBase = iP.Results.OutFileBase;

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
        [abfParams, abfData] = parse_abf(fileNames);

        % Extract the time, current and voltage vectors
        siMs = abfParams.siMs;
        tVec = abfData.tVec;
        iVecs = abfData.iVecs;
        vVecs = abfData.vVecs;
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

% Set output file base
if isempty(outFileBase)
    if iscell(fileBases)
        outFileBase = extract_common_prefix(fileBases);
    else
        outFileBase = fileBases;
    end
end

% Set output file base for titles
outFileBaseTitle = replace(outFileBase, '_', '\_');

% Create outFolder if not already exists
check_dir(outFolder);

% Count the number of sweeps
nSweeps = count_vectors(vVecs);

%% Parse the current injection protocol
% Identify the (common) current pulse endpoints
% TODO: Try using parse_pulse.m
[~, pulseEndPoints] = identify_CI_protocol(iVecs, siMs);

% Extract the current traces during the injection
iVecsInjected = extract_subvectors(iVecs, 'EndPoints', pulseEndPoints);

% Compute current injection values 
currentInjected = transpose(mean(iVecsInjected, 1));

%% Detect spikes
[spikesParams, spikesData] = detect_spikes_current_clamp(vVecs);

% Extract results
minAmpBefore = spikesParams.minAmpBefore;
minAmpAfter = spikesParams.minAmpAfter;
indSpikes = spikesData.indSpikes;

%% Compute spike frequencies during the injection
% Force as a column cell array
indSpikes = force_column_cell(indSpikes);

% Restrict to spikes occuring within the current pulse
indSpikesWithinPulse = extract_subvectors(indSpikes, 'Windows', pulseEndPoints);

% Compute spike frequency for each sweep
spikeFrequencyHz = compute_spike_frequency(indSpikesWithinPulse, siMs);

%% Return outputs as tables
% Force as a column cell array
[pulseEndPoints, timeVec, currentVec, voltageVec] = ...
    argfun(@force_column_cell, pulseEndPoints, tVec, iVecs, vVecs);

% Match row counts
[minAmpBefore, minAmpAfter, pulseEndPoints, ...
        timeVec, currentVec, voltageVec] = ...
    argfun(@(x) match_row_count(x, nSweeps), ...
            minAmpBefore, minAmpAfter, pulseEndPoints, ...
            timeVec, currentVec, voltageVec);

% Put together as tables
parsedParams = table(currentInjected, spikeFrequencyHz, ...
                    minAmpBefore, minAmpAfter);
parsedData = table(timeVec, currentVec, voltageVec, ...
                    pulseEndPoints, indSpikes, indSpikesWithinPulse);

%% Save tables
% Create full paths to files
sheetNameParams = fullfile(outFolder, [outFileBase, '_', ...
                            parsedParamsSuffix, '.csv']);
fileNameData = fullfile(outFolder, [outFileBase, '_', ...
                            parsedDataSuffix, '.mat']);

% Write the params into a spreadsheet
writetable(parsedParams, sheetNameParams);

% Save the data as a .mat file
save(fileNameData, 'parsedParams', 'parsedData', '-v7.3');

%% Plots

% Determine x and y limits
if plotRawData || plotSpikeDetection
    xLimits = compute_axis_limits(timeVec, 'x');
    yLimits = compute_axis_limits(voltageVec, 'y', 'Coverage', 90);
end

% Plot all sweeps together
if plotRawData
    % Plot all voltage traces
    figHandleVoltage = set_figure_properties('AlwaysNew', true);
    figNameVoltage = fullfile(outFolder, [outFileBase, '_raw_voltage_traces.png']);
    plot_traces(timeVec, voltageVec, 'PlotMode', 'overlapped', ...
                'XLabel', 'Time (ms)', 'YLabel', 'Membrane Potential (mV)', ...
                'FigHandle', figHandleVoltage, 'FigName', figNameVoltage);

    % Plot all current traces
    figHandleCurrent = set_figure_properties('AlwaysNew', true);
    figNameCurrent = fullfile(outFolder, [outFileBase, '_raw_current_traces.png']);
    plot_traces(timeVec, currentVec, 'PlotMode', 'overlapped', ...
                'XLabel', 'Time (ms)', 'YLabel', 'Current Injection (pA)', ...
                'FigHandle', figHandleCurrent, 'FigName', figNameCurrent);
end

% Plot each sweep individually with detected spikes and computed spike frequency    
if plotSpikeDetection
    % Create figure handles
    if plotSeparately
        % Create figure and axes handles
        [figHandlesSp, axHandlesSp] = ...
            arrayfun(@(x) create_subplots(1, 1, 'AlwaysNew', true), ...
                    transpose(1:nSweeps), 'UniformOutput', false);
    else
        % Create figure with subplots
        [figHandleSp, axHandlesSp] = ...
            create_subplots(nSweeps, 1, 'FigExpansion', [1, nSweeps / 3]);

        % Convert to cell array
        axHandlesSp = num2cell(axHandlesSp);
    end

    % Create figure names and titles
    if plotSeparately
        figNamesSp = create_labels_from_numbers(1:nSweeps, ...
                    'Prefix', fullfile(outFolder, [outFileBase, '_sweep']));

        % Create figure titles
        figTitlesSp = create_labels_from_numbers(1:nSweeps, ...
                    'Prefix', 'Spike Detection for Sweep #', ...
                    'Suffix', [' of ', outFileBaseTitle]);
    else
        figNameSp = fullfile(outFolder, [outFileBase, '_spike_detection.png']);
        figTitleSp = ['Spike Detection for ', outFileBaseTitle];
    end

    % Decide on x and y labels
    if plotSeparately
        xLabel = 'Time (ms)';
        yLabel = 'Membrane Potential (mV)';
    else
        xLabel = '';
        yLabel = '';
    end

    % Plot spike detection for each sweep separately
    cellfun(@(a, b, c, d, e, f, g) ...
                plot_spike_detection(a, b, c, d, e, f, g, ...
                                    xLimits, yLimits, xLabel, yLabel), ...
            timeVec, voltageVec, indSpikes, indSpikesWithinPulse, ...
            pulseEndPoints, num2cell(spikeFrequencyHz), axHandlesSp);

    % Save plots
    if plotSeparately
        for iSwp = 1:nSweeps
            figure(figHandlesSp{iSwp});
            title(figTitlesSp{iSwp});
            save_all_figtypes(figHandlesSp{iSwp}, figNamesSp{iSwp}, ...
                                {'png', 'epsc'})
        end
    else
        figure(figHandleSp);
%       suptitle(figTitleSp);
        save_all_figtypes(figHandleSp, figNameSp, {'png', 'epsc'})
    end
end

% Plot spike frequency over current injected (F-I plot)
if plotFI
    xLabel = 'Current Injected (pA)';
    yLabel = 'Spike Frequency (Hz)';
    figTitleFI = ['F-I plot for ', outFileBaseTitle];
    figHandleFI = set_figure_properties('AlwaysNew', true);
    figNameFI = fullfile(outFolder, [outFileBase, '_FI.png']);
    plot_tuning_curve(currentInjected, spikeFrequencyHz, ...
                'PLabel', xLabel, 'ReadoutLabel', yLabel, ...
                'FigTitle', figTitleFI,...
                'FigHandle', figHandleFI, 'FigName', figNameFI);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_spike_detection(tVec, vVec, indSpikes, indSpikesWithinPulse, ...
                                pulseEndPoints, spikeFrequencyHz, axHandle, ...
                                xLimits, yLimits, xLabel, yLabel)

% Go to appropriate axes
axes(axHandle);

% Hold on
hold on;

% Plot the voltage vector
plot(tVec, vVec, 'k');

% Plot spikes with red crosses
plot(tVec(indSpikes), vVec(indSpikes), 'rx');

% Plot spikes within pulse with green crosses
plot(tVec(indSpikesWithinPulse), vVec(indSpikesWithinPulse), 'gx');

% Set x and y limits
if ~isempty(xLimits)
    xlim(xLimits);
end
if ~isempty(yLimits)
    ylim(yLimits);
end

% Plot pulse boundaries
plot_window_boundaries(tVec(pulseEndPoints), 'BoundaryType', 'verticalLines');

% Annotate spike frequency
spikeFreqText = ['Spike Frequency: ', num2str(spikeFrequencyHz), ' Hz'];
text(0.2, 0.95, spikeFreqText, 'Units', 'normalized');

% Create labels
if ~isempty(xLabel)
    xlabel(xLabel);
end
if ~isempty(yLabel)
    ylabel(yLabel);
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

indSpikesWithinPulse = ...
    cellfun(@(x) x(x > pulseEndPoints(1) & x < pulseEndPoints(2)), ...
            indSpikes, 'UniformOutput', false);

h = figure(nSweeps + 1);
clf(h);
figName = [outFileBase, '_all'];
cm = colormap(lines);
for iSwp = 1:nSweeps
    plot(tVec, vVecs(:, iSwp), 'Color', cm(mod(iSwp, size(cm, 1)) + 1, :), ...
        'Displayname', ['Sweep #', num2str(iSwp)]);
    hold on;
end
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
legend('Location', 'northeast');
saveas(h, fullfile(outFolder, figName), 'png');

line([tVec(pulseEndPoints(1)), tVec(pulseEndPoints(1))], [y(1) y(2)]);
line([tVec(pulseEndPoints(2)), tVec(pulseEndPoints(2))], [y(1) y(2)]);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
