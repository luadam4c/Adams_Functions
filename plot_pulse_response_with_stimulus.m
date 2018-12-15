function h = plot_pulse_response_with_stimulus (tVec, respVec, stimVec, varargin)
%% Plots a pulse response with its stimulus
% Usage: h = plot_pulse_response_with_stimulus (tVec, respVec, stimVec, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - handle to figure
%                   specified as a figure handle
% Arguments:
%       tVec        - TODO: Description of reqarg1
%                   must be a numeric vector
%       respVec     - TODO: Description of reqarg1
%                   must be a numeric vector
%       stimVec     - TODO: Description of reqarg1
%                   must be a numeric vector
%       varargin    - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'OutFolder': the name of the directory for plots
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%                   - 'SaveFlag': whether to save the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'TODO': TODO
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/compute_and_plot_average_response.m

% File History:
% 2018-12-15 Moved from compute_and_plot_evoked_LFP.m
% 

%% Hard-coded parameters
colorAnnotations = 'r';

%% Default values for optional arguments
baseValueDefault = [];          % set later
peakValueDefault = [];          % set later
idxPeakDefault = [];            % set later
peakAmplitudeDefault = [];      % set later
labelsDefault = [];             % set later
outFolderDefault = pwd;         % save in present working directory by default
fileBaseDefault = 'Unnamed';
fileSuffixDefault = '_pulse_response';
responseNameDefault = 'Pulse Response';
saveFlagDefault = true;         % save the pulse train series by default
figTypesDefault = 'png';        % default figure type(s) for saving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'tVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'respVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'stimVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BaseValue', baseValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PeakValue', peakValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'IdxPeak', idxPeakDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'PeakAmplitude', peakAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Labels', labelsDefault, ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileSuffix', fileSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ResponseName', responseNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, tVec, respVec, stimVec, varargin{:});
baseValue = iP.Results.BaseValue;
peakValue = iP.Results.PeakValue;
idxPeak = iP.Results.IdxPeak;
peakAmplitude = iP.Results.PeakAmplitude;
labels = iP.Results.Labels;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
fileSuffix = iP.Results.FileSuffix;
responseName = iP.Results.ResponseName;
saveFlag = iP.Results.SaveFlag;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Preparation
% Check if needed output directories exist
check_dir(outFolder);

%% Do the job
% Open and clear figure
if saveFlag
    h = figure('Visible', 'off');
    figName = fullfile(outFolder, [fileBase, fileSuffix]);
    clf(h);
else
    figure;
end

% Compute the x axis limits
left = min(tVec);
right = max(tVec);
xLimits = [left, right];

% Find the time of the peak relative to the x limits
timeRange = right - left;
timePeakRel = (tVec(idxPeak) - left) / timeRange;

% Generate a subplot for the pulse response
%   Annotations:
%       red double arrow for peak amplitude
ax1 = subplot(3, 1, 1:2);
hold on;
plot(tVec, respVec);
xlim(xLimits);
ylimits = get(gca, 'YLim');
pos = get(gca, 'Position');
height = ylimits(2) - ylimits(1);
bottom = ylimits(1);
annotation('doublearrow', pos(1) + pos(3) * timePeakRel * ones(1, 2), ...
            pos(2) + pos(4) * ([baseValue, peakValue] - bottom) / height, ...
            'Color', colorAnnotations);
text(tVec(idxPeak) + 1, mean([baseValue, peakValue]), ...
        ['peak amp = ', num2str(peakAmplitude), ' (mV)']);
ylabel(labels{1});
title([responseName, ' for ', fileBase], 'Interpreter', 'none');

% Generate a subplot for the stimulation pulse
ax2 = subplot(3, 1, 3);
hold on;
plot(tVec, stimVec);
xlim(xLimits);    
ylabel(labels{2});
xlabel('Time (ms)');
title(['Stimulus for ', fileBase], 'Interpreter', 'none');

% Link the axes
linkaxes([ax1, ax2], 'x');

% Save and close figure
if saveFlag
    save_all_figtypes(h, figName, figTypes);
    close(h)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%