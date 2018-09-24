function [tVecLfp, vVecLfp, iVecStim, features] = compute_and_plot_evoked_LFP (fileName, varargin)
%% Computes and plots an evoked local field potential with its stimulus
% Usage: [tVecLfp, vVecLfp, iVecStim, features] = compute_and_plot_evoked_LFP (fileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       [tVecLfp, vVecLfp, iVecStim, features] = ...
%           compute_and_plot_evoked_LFP('20180914C_0001');
% Outputs:
%       tVecLfp     - time vector for evoked local field potential
%                   specified as a numeric column vector
%       vVecLfp     - voltage trace of evoked local field potential
%                   specified as a numeric column vector
%       iVecStim    - current trace of stimulation current pulse
%                   specified as a numeric column vector
%       features    - computed LFP features
%                       peakAmp
%                       peakSlope
% Arguments:    
%       fileName    - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       varargin    - 'OutFolder': the name of the directory for plots
%                   must be a string scalar or a character vector
%                   default == a subdirectory named by {fileName}_traces in pwd
%                   - 'PlotFlag': whether to plot the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveFlag': whether to save the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ChannelTypes': the channel types
%                   must be a cellstr with nChannels elements
%                       each being one of the following:
%                           'Voltage'
%                           'Current'
%                           'Conductance'
%                           'Undefined'
%                   default == detected with identify_channels()
%
% Requires:
%       cd/parse_abf.m
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%
% Used by:    
%       cd/plot_all_abfs_dir.m

% File History:
% 2018-09-17 Created by Adam Lu
% 2018-09-21 Considered the case when iVecs or vVecs is a cellarray
% 2018-09-21 Considered the case when iVecs or vVecs is 3-D
% 2018-09-23 Added computation of peak amplitude
% TODO: add timeUnits as a parameter with default 'ms'
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Undefined'};
baselineLengthMs = 5;           % baseline length in ms
colorAnnotations = 'r';

%% Default values for optional arguments
outFolderDefault = '';          % set later
plotFlagDefault = true;         % plot the evoked LFP with stim by default
saveFlagDefault = true;         % save the pulse train series by default
figTypesDefault = 'png';        % default figure type(s) for saving
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'fileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || iscellstr(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || iscellstr(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x));

% Read from the Input Parser
parse(iP, fileName, varargin{:});
outFolder = iP.Results.OutFolder;
plotFlag = iP.Results.PlotFlag;
saveFlag = iP.Results.SaveFlag;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

% Set (some) dependent argument defaults
[fileDir, fileBase, ~] = fileparts(fileName);
if isempty(outFolder)
    outFolder = fullfile(fileDir, strcat(fileBase, '_traces'));
end

%% Check if needed output directories exist
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end

%% Load data and prepare for plotting
% Load and parse the abf file
[abfParams, ~, tVec, vVecs, iVecs] = ...
    parse_abf(fileName, 'Verbose', false, ...
              'ChannelTypes', channelTypes, ...
              'ChannelUnits', channelUnits, ...
              'ChannelLabels', channelLabels);

% If vVecs or iVecs is a cellarray, use the first element
if iscell(vVecs)
    vVecs = vVecs{1};
end
if iscell(iVecs)
    iVecs = iVecs{1};
end

% If vVecs or iVecs is 3-D, use the first two dimensions 
if ndims(vVecs) > 2
    vVecs = squeeze(vVecs(:, :, 1));
end
if ndims(iVecs) > 2
    iVecs = squeeze(iVecs(:, :, 1));
end

% Extract the parsed parameters
channelLabels = abfParams.channelLabels;
nSweeps = abfParams.nSweeps;
siMs = abfParams.siMs;

% Compute the baseline length in samples
baselineLengthSamples = floor(baselineLengthMs / siMs);

%% Average the voltage responses
% Identify the current pulse response endpoints
idxCprStarts = zeros(nSweeps, 1);
idxCprEnds = zeros(nSweeps, 1);
idxCpStarts = zeros(nSweeps, 1);
idxCpEnds = zeros(nSweeps, 1);
parfor iSwp = 1:nSweeps
    % Extract the voltage and current vectors for this sweep
    vVecCpr = vVecs(:, iSwp);
    iVecCpr = iVecs(:, iSwp);

    % Identify the current pulse and current pulse response endpoints
    [idxCprStarts(iSwp), idxCprEnds(iSwp), ~, ...
        idxCpStarts(iSwp), idxCpEnds(iSwp)] = ...
        find_pulse_response_endpoints(vVecCpr, siMs, ...
                                        'IvecCpr', iVecCpr, ...
                                        'BaselineLengthMs', baselineLengthMs);
end

% Compute the number of samples in each current pulse response
nSamplesEachCpr = idxCprEnds - idxCprStarts + 1;

% Determine number of samples in the shortest current pulse response
nSamplesCpr = min(nSamplesEachCpr);

% Extract the time vector using the starting index from the first sweep
idxCprStartFirst = idxCprStarts(1);
idxCprStartEnd = (idxCprStartFirst - 1) + nSamplesCpr;
tVecLfp = tVec(idxCprStartFirst:idxCprStartEnd);

% Place the current pulses and current pulse responses in the same matrix
iVecCprs = zeros(nSamplesCpr, nSweeps);
vVecCprs = zeros(nSamplesCpr, nSweeps);
parfor iSwp = 1:nSweeps
    % Get the starting index of the current pulse response for this sweep
    idxCprStart = idxCprStarts(iSwp);

    % Get the ending index of the current pulse response for this sweep
    idxCprEnd = (idxCprStart - 1) + nSamplesCpr;

    % Extract the current pulse and current pulse responses
    iVecCprs(:, iSwp) = iVecs(idxCprStart:idxCprEnd, iSwp);
    vVecCprs(:, iSwp) = vVecs(idxCprStart:idxCprEnd, iSwp);
end

% Average the current pulses to get the stimulation pulse
iVecStim = mean(iVecCprs, 2);

% Average the current pulse responses to get the evoked local field potential
vVecLfp = mean(vVecCprs, 2);

%% Extract the amplitude of the evoked local field potential
% Compute the baseline voltage value of the averaged trace
baseVal = mean(vVecLfp(1:baselineLengthSamples));

% Assume the current pulses all end at the same index
idxCpEnd = idxCpEnds(1) - idxCprStartFirst + 1;

% Compute the peak voltage value of the averaged trace
[peakVal, temp1] = max(vVecLfp((idxCpEnd + 1):end));
idxPeak = temp1 + idxCpEnd;

% Compute the relative peak amplitude
peakAmp = peakVal - baseVal;

%% Extract the slope of the evoked local field potential
% TODO

%% Plot the evoked local field potential with the stimulation pulse
if plotFlag
    % Open and clear figure
    if saveFlag
        h = figure('Visible', 'off');
        figName = fullfile(outFolder, [fileBase, '_LFP']);
        clf(h);
    else
        figure;
    end

    % Compute the x axis limits
    left = min(tVecLfp);
    right = max(tVecLfp);
    xlimits = [left, right];
    
    % Find the time of the peak relative to the x limits
    width = right - left;
    timePeakRel = (tVecLfp(idxPeak) - left) / width;

    % Generate a subplot for the evoked local field potential
    %   Annotations:
    %       red double arrow for peak amplitude
    ax1 = subplot(3, 1, 1:2);
    hold on;
    plot(tVecLfp, vVecLfp);
    xlim(xlimits);
    ylimits = get(gca, 'YLim');
    pos = get(gca, 'Position');
    height = ylimits(2) - ylimits(1);
    bottom = ylimits(1);
    annotation('doublearrow', pos(1) + pos(3) * timePeakRel * ones(1, 2), ...
                pos(2) + pos(4) * ([baseVal, peakVal] - bottom) / height, ...
                'Color', colorAnnotations);
    text(tVecLfp(idxPeak) + 1, mean([baseVal, peakVal]), ...
            ['peak amp = ', num2str(peakAmp), ' (mV)']);
    ylabel(channelLabels{1});
    title(['Evoked potential for ', fileBase], 'Interpreter', 'none');

    % Generate a subplot for the stimulation pulse
    ax2 = subplot(3, 1, 3);
    hold on;
    plot(tVecLfp, iVecStim);
    xlim(xlimits);    
    ylabel(channelLabels{2});
    xlabel('Time (ms)');
    title(['Stimulus for ', fileBase], 'Interpreter', 'none');

    % Link the axes
    linkaxes([ax1, ax2], 'x');

    % Save and close figure
    if saveFlag
        save_all_figtypes(h, figName, figTypes);
        close(h)
    end
end

% Output features in the features structure
features.peakAmp = peakAmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

idxCpStarts = zeros(nSweeps, 1);
idxCpEnds = zeros(nSweeps, 1);

[idxCprStart(iSwp), idxCprEnds(iSwp), ~, ...
    idxCpStarts(iSwp), idxCpEnds(iSwp)] = ...
    find_pulse_response_endpoints(vVecs(:, iSwp), siMs, ...
                                    'IvecCpr', iVecs(:, iSwp));

%}