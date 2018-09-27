function compute_and_plot_concatenated_trace(abfParamsAllStruct, dataReorderedAll, varargin)
%% Computes and plots concatenated traces from parsed ABF file results
% Usage: compute_and_plot_concatenated_trace(abfParamsAllStruct, dataReorderedAll, varargin)
%
%
% Used by:
%       cd/plot_all_abfs.m

% File History:
% 2018-09-21 - Moved from plot_all_abfs_dir.m 
% TODO: make the width dependent on timeLength

paperPosition = [0, 0, 80, 5];

sourceDirectoryDefault = pwd;
outFolderDefault = '';
outputLabelDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SourceDirectory', sourceDirectoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutputLabel', outputLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
sourceDirectory = iP.Results.SourceDirectory;
outFolder = iP.Results.OutFolder;
outputLabel = iP.Results.OutputLabel;

%% Preparation
% Decide on the outFolder
if isempty(outFolder)
    outFolder = sourceDirectory;
end

% Decide on the output label
[~, directoryName, ~] = fileparts(sourceDirectory);
if isempty(outputLabel)
    outputLabel = directoryName;
end
figName = fullfile(outFolder, [directoryName, '_combined']);

% Count the number of files
nFiles1 = numel(dataReorderedAll);
nFiles2 = numel(abfParamsAllStruct);
if nFiles1 ~= nFiles2
    fprintf(['abfParamsAllStruct and dataReorderedAll does not ', ...
                'have the same number of elements!\n'])
    fprintf('Sweeps cannot be concatenated!!\n\n');
    return;
else
    nFiles = nFiles1;
end

% TODO: make sure abfParamsAllStruct and dataReorderedAll have the same numel

% Check if the number of channels are all the same
nChannelsAll = [abfParamsAllStruct.nChannels];
if numel(unique(nChannelsAll)) ~= 1
    fprintf('The number of channels are not all the same!\n');
    fprintf('Sweeps cannot be concatenated!!\n\n');
    return;
end

% Check if the channel types are all the same
channelTypesAll = {abfParamsAllStruct.channelTypes};
if ~isequal(channelTypesAll{:})
    fprintf('The channel types are not all the same!\n');
    fprintf('Sweeps cannot be concatenated!!\n\n');
    return;
end

%% Concatenate sweeps within each file
nSamplesAll = zeros(nFiles, 1);
channelLabelsAll = cell(nFiles, 1);
tVecsAll = cell(nFiles, 1);
dataConcatenatedAll = cell(nFiles, 1);
parfor iFile = 1:nFiles
    % Extract from arrays
    abfParams = abfParamsAllStruct(iFile);
    dataReordered = dataReorderedAll{iFile};

    % Extract from structure
    nDimensions = abfParams.nDimensions;
    nSamples = abfParams.nSamples;
    nSweeps = abfParams.nSweeps;
    siSeconds = abfParams.siSeconds;
    channelLabels = abfParams.channelLabels;

    % Concatenate all sweeps if dimensions is greater than 2
    if nDimensions > 2
        % Get the number of samples, sweeps and channels
        dimsNew = size(dataReordered);

        % Compute the new number of samples
        nSamplesNew = nSamples * nSweeps; 
        dimsNew(1) = nSamplesNew;

        % Compute the new number of sweeps
        dimsNew(2) = 1;

        % Reshape the matrix
        dataConcatenated = squeeze(reshape(dataReordered, dimsNew));
    else
        % There is no need to concatenate
        dataConcatenated = dataReordered;
        nSamplesNew = nSamples;
    end

    % Construct a time vector in seconds
    tVec = (1:nSamplesNew)' * siSeconds;

    % Store in arrays
    nSamplesAll(iFile) = nSamplesNew;
    channelLabelsAll{iFile} = channelLabels;
    tVecsAll{iFile} = tVec;
    dataConcatenatedAll{iFile} = dataConcatenated;
end

%% Concatenate traces from all files
dataCombined = vertcat(dataConcatenatedAll{:});

%% Generate a combined time vector
% Compute the total number of samples
nSamplesCombined = sum(nSamplesAll);

tLast = 0;
% TODO:tVecCombined = zeros(nSamplesCombined);
tVecCombined = [];
for iFile = 1:nFiles
    % Get the current time vector
    tVecNow = tVecsAll{iFile} + tLast;

    % Combine with the previous vectors
    tVecCombined = [tVecCombined; tVecNow];

    % Update the last time
    tLast = tVecCombined(end);
end

%% Prepare for plotting
% Use the first channel labels
channelLabels = abfParamsAllStruct(1).channelLabels;
nChannels = abfParamsAllStruct(1).nChannels;

% Get the xlimits
xlimits = [min(tVecCombined), max(tVecCombined)];

%% Plot the combined trace
% Create and clear the figure
h = figure('Visible', 'off', 'PaperPosition', paperPosition);
clf(h)

% Plot each channel in a subplot separately
for iChannel = 1:nChannels
    % Create the subplot
    ax(iChannel) = subplot(nChannels, 1, iChannel);

    % Plot the vector
    plot(tVecCombined, dataCombined(:, iChannel));

    % Generate a y label
    ylabel(channelLabels{iChannel});

    % Define x limits
    if iChannel == nChannels
        xlabel('Time (seconds)');
    end

    % Define x limits
    xlim(xlimits);
end

% Link the x axes
linkaxes(ax, 'x');

% Don't use subscripts
outputLabel = strrep(outputLabel, '_', '\_');

% Plot a title
suptitle(['Combined traces for ', outputLabel]);

% Save and close the figure
print(h, figName, '-djpeg');
close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
