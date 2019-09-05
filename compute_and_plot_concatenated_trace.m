function compute_and_plot_concatenated_trace(abfParamsAllStruct, dataReorderedAll, varargin)
%% Computes and plots concatenated traces from parsed ABF file results
% Usage: compute_and_plot_concatenated_trace(abfParamsAllStruct, dataReorderedAll, varargin)
% Explanation:
%       Concatenates all sweeps of each channel in a directory of .abf files
% Outputs:
%       
% Arguments:
%       abfParamsAllStruct      - abf params of all files to be concatenated
%                               specified as a struct
%       dataReorderedAll        - 
%       varargin                
%                               - 'XLimits': limits of x axis
%                                           suppress by setting value to 'suppress'
%                               must be 'suppress' or a 2-element increasing numeric vector
%                               default == [min(tVecCombined), max(tVecCombined)]
%                               - 'XLabel': label for the time axis, 
%                                           suppress by setting value to 'suppress'
%                               must be a string scalar or a character vector 
%                                   or a cell array of strings or character vectors
%                               default == 'Time (seconds)'
%                               - 'YLabel': label(s) for the y axis, 
%                                           suppress by setting value to 'suppress'
%                               must be a string scalar or a character vector
%                               default == abfParamsAllStruct(1).channelLabels
%                               - Any other parameter-value pair for 
%                               the plot_traces() function
% Requires:
%
% Used by:
%       cd/plot_all_abfs.m

% File History:
% 2018-09-21 - Moved from plot_all_abfs_dir.m 
% 2019-04-15 BT - Implemented plot_traces.m
% TODO: make the width dependent on timeLength

paperPosition = [0, 0, 80, 5];

sourceDirectoryDefault = pwd;
outFolderDefault = '';
outputLabelDefault = '';
xLimitsDefault = [];
xLabelDefault = 'Time (minutes)';
yLabelDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SourceDirectory', sourceDirectoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutputLabel', outputLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || iscell(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, varargin{:});
sourceDirectory = iP.Results.SourceDirectory;
outFolder = iP.Results.OutFolder;
outputLabel = iP.Results.OutputLabel;
xLimits = iP.Results.XLimits;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;

% Keep unmatched arguments for the plot_traces() function
otherArguments = struct2arglist(iP.Unmatched);

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

% Decide on the y label
if isempty(yLabel)
    yLabel = abfParamsAllStruct(1).channelLabels;
end

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

% Check if the number of channels are all the same
nChannelsAll = [abfParamsAllStruct.nChannels];
if numel(unique(nChannelsAll)) ~= 1
    fprintf('The number of channels are not all the same!\n');
    fprintf('Sweeps cannot be concatenated!!\n\n');
    return;
end

% Check if the channel types are all the same
channelTypesAll = {abfParamsAllStruct.channelTypes};
if numel(channelTypesAll) > 2 && ~isequal(channelTypesAll{:})
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

% Convert tVecCombined to minutes
tVecCombined = tVecCombined / 60;

% Get the xlimits
if isempty(xLimits)
    xLimits = [min(tVecCombined), max(tVecCombined)];
end

%% Plot the combined trace
% Create and clear the figure
h = figure('Visible', 'off', 'PaperPosition', paperPosition);
clf(h)

outputLabel = strrep(outputLabel, '_', '\_');

h = plot_traces(tVecCombined, dataCombined, 'FigHandle', h, 'LinkAxesOption', 'x', ...
                'FigTitle', ['Combined traces for ', outputLabel], 'PlotMode', 'parallel', ...
                'XLimits', xLimits, 'XLabel', xLabel, 'YLabel', yLabel, ...
                'LegendLocation', 'suppress', ...
                otherArguments{:});

% Save and close the figure
print(h, figName, '-djpeg');
close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:
% Use the first channel labels
channelLabels = abfParamsAllStruct(1).channelLabels;
nChannels = abfParamsAllStruct(1).nChannels;

% % Plot each channel in a subplot separately
% for iChannel = 1:nChannels
%     % % Create the subplot
%     % ax(iChannel) = subplot(nChannels, 1, iChannel);

%     % Plot the vector
%     % TODO: Use plot_traces.m instead
%     % Identify if arguments to plot_traces were passed
%     xlim_idx = find(strcmp(otherArguments, 'XLimits'));
%     xlabel_idx = find(strcmp(otherArguments, 'XLabel'));
%     ylabel_idx = find(strcmp(otherArguments, 'YLabel'));
%     % If no arguments passed, use default
%     if isempty(xlim_idx)
%         otherArguments{end+1} = 'XLimits'; otherArguments{end+1} = xlimits;
%     end
%     if isempty(xlabel_idx)
%         if iChannel == nChannels
%             otherArguments{end+1} = 'XLabel'; otherArguments{end+1} = 'Time (seconds)';
%         end
%     end
%     if isempty(ylabel_idx)
%         otherArguments{end+1} = 'YLabel'; otherArguments{end+1} = channelLabels{iChannel};
%     end
%     [h, ax(iChannel)] = plot_traces(tVecCombined, dataCombined(:, iChannel), 'FigHandle', h, 'YLimits', [-0.1 0.1], otherArguments{:});
%     %{
%     plot(tVecCombined, dataCombined(:, iChannel));

%     % Generate a y label
%     ylabel(channelLabels{iChannel});

%     % Define x limits
%     if iChannel == nChannels
%         xlabel('Time (seconds)');
%     end

%     % Define x limits
%     xlim(xlimits);
%     %}
% end

% % Link the x axes
% linkaxes(ax, 'x');

% % Don't use subscripts
% outputLabel = strrep(outputLabel, '_', '\_');

% % Plot a title
% suptitle(['Combined traces for ', outputLabel]);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%