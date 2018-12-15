function [tVecLfp, vVecLfp, iVecStim, features, h] = ...
                compute_and_plot_evoked_LFP (fileName, varargin)
%% Computes and plots an evoked local field potential with its stimulus
% Usage: [tVecLfp, vVecLfp, iVecStim, features, h] = ...
%               compute_and_plot_evoked_LFP (fileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       [tVecLfp, vVecLfp, iVecStim, features, h] = ...
%           compute_and_plot_evoked_LFP('20180914C_0001');
% Outputs:
%       tVecLfp     - time vector for evoked local field potential
%                   specified as a numeric column vector
%       vVecLfp     - voltage trace of evoked local field potential
%                   specified as a numeric column vector
%       iVecStim    - current trace of stimulation current pulse
%                   specified as a numeric column vector
%       features    - computed LFP features, a table
%                       peakAmplitude
%                       peakSlope
%       h           - handle to figure
%                   specified as a figure handle
% Arguments:    
%       fileName    - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       varargin    - 'OutFolder': the name of the directory for plots
%                   must be a string scalar or a character vector
%                   default == a subdirectory named by outFolderName
%                               in the same directory as fileName
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
%                           'Other'
%                   default == detected with identify_channels()
%                   - 'ParsedParams': parsed parameters returned by parse_abf.m
%                   must be a scalar structure
%                   default == what the file provides
%                   - 'ParsedData': parsed data returned by parse_abf.m
%                   must be a scalar structure
%                   default == what the file provides
%
% Requires:
%       cd/compute_average_pulse_response.m
%       cd/isfigtype.m
%       cd/plot_pulse_response_with_stimulus.m
%       cd/save_all_figtypes.m
%
% Used by:    
%       cd/plot_all_abfs.m

% File History:
% 2018-09-17 - Created by Adam Lu
% 2018-09-21 - Considered the case when iVecs or vVecs is a cellarray
% 2018-09-21 - Considered the case when iVecs or vVecs is 3-D
% 2018-09-23 - Added computation of peak amplitude
% 2018-09-30 - Changed default outFolder to 'LFPs' under fileDir
% 2018-10-03 - Updated usage of parse_abf.m
% 2018-10-03 - Added ParsedData, ParsedParams as optional arguments
% 2018-10-09 - Updated usage of find_pulse_response_endpoints.m
% 2018-12-15 - Updated usage of find_pulse_response_endpoints.m
% 2018-12-15 - Now uses argfun.m, compute_average_trace.m, 
%               extract_subvectors.m, force_column_cell.m
% 2018-12-15 Moved code to compute_and_plot_average_pulse_response.m
% TODO: add timeUnits as a parameter with default 'ms'
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

% For computing
responseType = 'Voltage';
lowPassFrequency = 1000;        % lowpass filter frequency in Hz
baselineLengthMs = 5;           % baseline length in ms
responseLengthMs = 20;          % response length in ms

% For plotting
outFolderName = 'LFPs';
fileSuffix = '_LFP';
responseName = 'Evoked potential';

%% Default values for optional arguments
outFolderDefault = '';          % set later
plotFlagDefault = true;         % plot the evoked LFP with stim by default
saveFlagDefault = true;         % save the pulse train series by default
figTypesDefault = 'png';        % default figure type(s) for saving
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later
parsedParamsDefault = [];       % set later
parsedDataDefault = [];         % set later

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
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ParsedParams', parsedParamsDefault, ...
    @(x) isempty(x) || isstruct(x));
addParameter(iP, 'ParsedData', parsedDataDefault, ...
    @(x) isempty(x) || isstruct(x));

% Read from the Input Parser
parse(iP, fileName, varargin{:});
outFolder = iP.Results.OutFolder;
plotFlag = iP.Results.PlotFlag;
saveFlag = iP.Results.SaveFlag;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;
parsedParams = iP.Results.ParsedParams;
parsedData = iP.Results.ParsedData;

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

%% Rename outputs
[tVecLfp, vVecLfp, iVecStim, features, h] = ...
    compute_and_plot_average_pulse_response (fileName, responseType, ...
        'LowPassFrequency', lowPassFrequency, ...
        'BaselineLengthMs', baselineLengthMs, ...
        'ResponseLengthMs', responseLengthMs, ...
        'OutFolder', outFolder, 'OutFolderName', outFolderName, ...
        'FileSuffix', fileSuffix, 'ResponseName', responseName, ...
        'PlotFlag', plotFlag, 'SaveFlag', saveFlag, 'FigTypes', figTypes, ...
        'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
        'ChannelLabels', channelLabels, ...
        'ParsedParams', parsedParams, 'ParsedData', parsedData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

idxCpStarts = zeros(nSweeps, 1);
idxCpEnds = zeros(nSweeps, 1);

[idxCprStart(iSwp), idxCprEnds(iSwp), ~, ...
    idxCpStarts(iSwp), idxCpEnds(iSwp)] = ...
    find_pulse_response_endpoints(vVecs(:, iSwp), siMs, ...
                                    'IvecCpr', iVecs(:, iSwp));

outFolder = fullfile(fileDir, strcat(fileBase, '_traces'));

[parsedParams, ~, tVec, vVecs, iVecs] = ...
    parse_abf(fileName, 'Verbose', false, ...
              'ChannelTypes', channelTypes, ...
              'ChannelUnits', channelUnits, ...
              'ChannelLabels', channelLabels);

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
                                        'PulseVectors', iVecCpr, ...
                                        'BaselineLengthMs', baselineLengthMs, ...
                                        'ResponseLengthMs', responseLengthMs);
end

% Compute the number of samples in each current pulse response
nSamplesEachCpr = idxCprEnds - idxCprStarts + 1;

% Determine number of samples in the shortest current pulse response
nSamplesCpr = min(nSamplesEachCpr);

% Extract the time vector using the starting index from the first sweep
idxCprStartFirst = idxCprStarts(1);
idxCprEndFirst = (idxCprStartFirst - 1) + nSamplesCpr;
tVecLfp = tVec(idxCprStartFirst:idxCprEndFirst);

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

idxCpEnd = idxCpEnds(1) - idxCprStartFirst + 1;

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

nSweeps = parsedParams.nSweeps;

% Compute the baseline voltage value of the averaged trace
baseValue = mean(vVecLfp(1:baselineLengthSamples));

% Get the index of pulse end for the restricted trace
%   Note: Assume the pulses all end at the same index
idxResponseEnd = idxCpEnds(1) - idxCprStartAveraged + 1;

% Compute the peak voltage value of the averaged pulse response trace
%   after the pulse ends
[peakValue, temp1] = max(vVecLfp((idxResponseEnd + 1):end));
idxPeak = temp1 + idxResponseEnd;

% Compute the relative peak amplitude
peakAmplitude = peakValue - baseValue;

% Output features in the features structure
features.peakAmplitude = peakAmplitude;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%