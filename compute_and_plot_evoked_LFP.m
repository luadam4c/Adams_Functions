function [tVecLfp, vVecLfp, iVecStim] = compute_and_plot_evoked_LFP (fileName, varargin)
%% Plots an evoked local field potential
% Usage: [tVecLfp, vVecLfp, iVecStim] = compute_and_plot_evoked_LFP (fileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       [tVecLfp, vVecLfp, iVecStim] = ...
%           compute_and_plot_evoked_LFP('20180914C_0001');
% Outputs:
%       tVecLfp     - time vector for evoked local field potential
%                   specified as a numeric column vector
%       vVecLfp     - voltage trace of evoked local field potential
%                   specified as a numeric column vector
%       iVecStim      - current trace of stimulation current pulse
%                   specified as a numeric column vector
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
% 

%% Hard-coded parameters

%% Default values for optional arguments
outFolderDefault = '';          % set later
plotFlagDefault = true;         % plot the evoked LFP with stim by default
saveFlagDefault = true;         % save the pulse train series by default
figTypesDefault = 'png';        % default figure type(s) for saving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
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

% Read from the Input Parser
parse(iP, fileName, varargin{:});
outFolder = iP.Results.OutFolder;
plotFlag = iP.Results.PlotFlag;
saveFlag = iP.Results.SaveFlag;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

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
[abfParams, ~, tVec, vVecs, iVecs] = parse_abf(fileName, 'Verbose', false);

% Extract the parsed parameters
channelLabels = abfParams.channelLabels;
nSweeps = abfParams.nSweeps;
siMs = abfParams.siMs;

%% Average the voltage responses
% Identify the current pulse response endpoints
idxCprStarts = zeros(nSweeps, 1);
idxCprEnds = zeros(nSweeps, 1);
parfor iSwp = 1:nSweeps
    % Identify the current pulse and current pulse response endpoints
    [idxCprStarts(iSwp), idxCprEnds(iSwp)] = ...
        find_pulse_response_endpoints(vVecs(:, iSwp), siMs, ...
                                        'IvecCpr', iVecs(:, iSwp));
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
%   and extract the time vector from the first response
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

    % Generate a subplot for the evoked local field potential
    ax1 = subplot(3, 1, 1:2);
    plot(tVecLfp, vVecLfp);
    ylabel(channelLabels{1});       % TODO: Make more robust
    title(['Evoked potential for ', fileBase], 'Interpreter', 'none');

    % Generate a subplot for the stimulation pulse
    ax2 = subplot(3, 1, 3);
    plot(tVecLfp, iVecStim);
    ylabel(channelLabels{2});       % TODO: Make more robust
    xlabel('Time (ms)');
    title(['Stimlus for ', fileBase], 'Interpreter', 'none');

    % Link the axes
    linkaxes([ax1, ax2], 'x');

    % Adjust the x limits
    xlim([min(tVecLfp), max(tVecLfp)]);

    % Save and close figure
    if saveFlag
        save_all_figtypes(h, figName, figTypes);
        close(h)
    end
end

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