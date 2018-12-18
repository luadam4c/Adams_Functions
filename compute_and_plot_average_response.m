function [tVecAvg, respAvg, stimAvg, featuresAvg, h] = ...
            compute_and_plot_average_response (fileName, responseType, varargin)
%% Computes and plots an average pulse response with its stimulus
% Usage: [tVecAvg, respAvg, stimAvg, featuresAvg, h] = ...
%           compute_and_plot_average_response (fileName, responseType, varargin)
% Explanation:
%       TODO
% Example(s):
%       [tVecAvg, respAvg, stimAvg, featuresAvg, h] = ...
%               compute_and_plot_average_response('20180914C_0001', 'Voltage');
% Outputs:
%       tVecAvg     - time vector for average response
%                   specified as a numeric column vector
%       respAvg     - average pulse response vector
%                   specified as a numeric column vector
%       stimAvg     - average stimulation pulse vector
%                   specified as a numeric column vector
%       featuresAvg - computed features of the average pulse response, 
%                       returned by compute_average_pulse_response.m
%                   specified as a table of 1 row
%       h           - handle to figure
%                   specified as a figure handle
% Arguments:    
%       fileName    - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       responseType- the channel type for the pulse response
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'Voltage'       - voltage
%                       'Current'       - current
%                       'Conductance'   - conductance
%                       'Other'         - other un-identified types
%       varargin    - 'LowPassFrequency': frequency of lowpass filter in Hz
%                   must be a nonnegative scalar
%                   default = []
%                   - 'ResponseLengthMs': length of the pulse response
%                                           after pulse endpoint in ms
%                   must be a nonnegative scalar
%                   default = 20 ms
%                   - 'BaselineLengthMs': length of the pulse response
%                                           before pulse start in ms
%                   must be a nonnegative scalar
%                   default = 5 ms
%                   - 'MinPeakDelayMs': minimum peak delay (ms)
%                               after the end of the pulse
%                   must be a positive scalar
%                   default == 0 ms
%                   - 'OutFolder': the name of the directory for plots
%                   must be a string scalar or a character vector
%                   default == a subdirectory named by outFolderName
%                               in the same directory as fileName
%                   - 'OutFolderName': the name of the default output directory
%                   must be a string scalar or a character vector
%                   default == 'Pulse_Response'
%                   - 'FileSuffix': file suffix
%                   must be a string scalar or a character vector
%                   default == '_pulse_response'
%                   - 'ResponseName': name of current pulse response
%                   must be a string scalar or a character vector
%                   default == 'Pulse Response'
%                   - 'PlotFlag': whether to plot the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SavePlotsFlag': whether to save plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveTablesFlag': whether to save tables
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
%       cd/plot_protocols.m

% File History:
% 2018-12-15 - Moved from compute_and_plot_evoked_LFP.m
% 2018-12-15 - Now used directly in plot_protocols.m
% 2018-12-15 - Now uses the file name as the row name for the feature table

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

%% Default values for optional arguments
lowPassFrequencyDefault = [];   % do not lowpass filter by default
responseLengthMsDefault = 20;   % a response of 20 ms by default
baselineLengthMsDefault = 5;    % a baseline of 5 ms by default
minPeakDelayMsDefault = 0;      % no minimum peak delay by default
outFolderDefault = '';          % set later
outFolderNameDefault = 'Pulse_Response';
fileSuffixDefault = '_pulse_response';
responseNameDefault = 'Pulse Response';
plotFlagDefault = true;         % plot the evoked LFP with stim by default
savePlotsFlagDefault = true;    % save all plots by default
saveTablesFlagDefault = true;   % save all tables by default
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
addRequired(iP, 'responseType', ...
    @(x) any(validatestring(x, validChannelTypes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LowPassFrequency', lowPassFrequencyDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'ResponseLengthMs', responseLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'BaselineLengthMs', baselineLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolderName', outFolderNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileSuffix', fileSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ResponseName', responseNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SavePlotsFlag', savePlotsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveTablesFlag', saveTablesFlagDefault, ...
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
parse(iP, fileName, responseType, varargin{:});
lowPassFrequency = iP.Results.LowPassFrequency;
responseLengthMs = iP.Results.ResponseLengthMs;
baselineLengthMs = iP.Results.BaselineLengthMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
outFolder = iP.Results.OutFolder;
outFolderName = iP.Results.OutFolderName;
fileSuffix = iP.Results.FileSuffix;
responseName = iP.Results.ResponseName;
plotFlag = iP.Results.PlotFlag;
savePlotsFlag = iP.Results.SavePlotsFlag;
saveTablesFlag = iP.Results.SaveTablesFlag;
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

% Set (some) dependent argument defaults
[fileDir, fileBase, ~] = fileparts(fileName);
if isempty(outFolder)
    outFolder = fullfile(fileDir, outFolderName);
end

%% Average the pulse and pulse responses
[tVecAvg, respAvg, stimAvg, featuresAvg] = ...
    compute_average_pulse_response(fileName, responseType, ...
        'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
        'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
        'ChannelLabels', channelLabels, ...
        'LowPassFrequency', lowPassFrequency, ...
        'BaselineLengthMs', baselineLengthMs, ...
        'ResponseLengthMs', responseLengthMs, ...
        'MinPeakDelayMs', minPeakDelayMs, ...
        'SaveFlag', saveTablesFlag);

%% Plot the evoked local field potential with the stimulation pulse
if plotFlag
    % Save in a single params structure
    params = table2struct(featuresAvg);
    params.OutFolder = outFolder;
    params.SaveFlag = savePlotsFlag;
    params.FigTypes = figTypes;
    params.FileBase = fileBase;
    params.FileSuffix = fileSuffix;
    params.ResponseName = responseName;

    % Plot the pulse response with the stimulation pulse
    h = plot_pulse_response_with_stimulus(tVecAvg, respAvg, stimAvg, params);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%