function [abfParamsAllStruct, dataAll, tVecAll, vVecsAll, iVecsAll, ...
            gVecsAll, dataReorderedAll, abfParamsAllCell] = ...
                parse_all_abfs (varargin)
%% Parses all abf files in the directory
% Usage: [abfParamsAllStruct, dataAll, tVecAll, vVecsAll, iVecsAll, ...
%           gVecsAll, dataReorderedAll, abfParamsAllCell] = ...
%               parse_all_abfs (varargin)
% Explanation:
%       This function calls parse_abf.m with 'IdentifyProtocols' == true
%           for all the .abf files in the provided directory (default pwd)
% Example(s):
%       [abfParams, data, tVec, vVecs, iVecs, gVecs] = ...
%           parse_all_abfs;
% Outputs:
%       structure array of abfParams (see parse_abf.m)
%       cell arrays of other outputs of parse_abf.m
%
% Arguments:
%       varargin    - 'Directory': the name of the directory containing 
%                                   the abf files, e.g. '20161216'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileNames': names of .abf files to detect
%                   must be a cell array of character arrays or strings
%                   default == detect from pwd
%                   - 'Verbose': whether to output parsed results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'UseOriginal': whether to use original 
%                           channel labels and units over identify_channels()
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ExpMode': experiment mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'EEG'   - EEG data; x axis in seconds; y-axis in uV
%                       'patch' - patch data; x axis in ms; y-axis in mV
%                   default == 'EEG' for 2d data 'patch' for 3d data
%                   - 'TimeUnits': units for time
%                   must be a string scalar or a character vector
%                   default == 's' for 2-data data and 'ms' for 3-data data
%                   - 'ChannelTypes': the channel types
%                   must be a cellstr with nChannels elements
%                       each being one of the following:
%                           'Voltage'
%                           'Current'
%                           'Conductance'
%                           'Undefined'
%                   default == detected with identify_channels()
%                   - 'ChannelUnits': the channel units
%                   must be a cellstr with nChannels elements
%                   default == detected with identify_channels()
%                   - 'ChannelLabels': the channel labels
%                   must be a cellstr with nChannels elements
%                   default == detected with identify_channels()
%                   - 'IdentifyProtocols': whether to identify protocols
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/parse_abf.m
%       /home/Matlab/Downloaded_Functions/dirr.m
%
% Used by:    
%       cd/plot_all_abfs.m

% File History:
% 2018-09-27 Pulled code from plot_all_abfs.m
% 

%% Hard-coded parameters
validExpModes = {'EEG', 'patch', ''};
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Undefined'};

%% Default values for optional arguments
directoryDefault = pwd;             % look for .abf files in 
                                    %   the present working directory by default
fileNamesDefault = {};              % detect from pwd by default
verboseDefault = false;             % print to standard output by default
useOriginalDefault = false;         % use identify_channels.m instead
                                    % of the original channel labels by default
expModeDefault = '';                % set later
timeUnitsDefault = '';              % set later
channelTypesDefault = {};           % set later
channelUnitsDefault = {};           % set later
channelLabelsDefault = {};          % set later
identifyProtocolsDefault = false;   % don't identify protocols by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'UseOriginal', useOriginalDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ExpMode', expModeDefault, ...
    @(x) isempty(x) || any(validatestring(x, validExpModes)));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'IdentifyProtocols', identifyProtocolsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
verbose = iP.Results.Verbose;
useOriginal = iP.Results.UseOriginal;
expMode = validatestring(iP.Results.ExpMode, validExpModes);
timeUnits = iP.Results.TimeUnits;
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;
identifyProtocols = iP.Results.IdentifyProtocols;

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

%% Get file names
% Decide on the files to use
if isempty(fileNames)
    % Find all .abf files in the directory
    [~, ~, fileNames] = dirr(directory, '.abf', 'name');
    if isempty(fileNames)
        fprintf('No abf files in current directory!\n');
        fprintf('Type ''help plot_all_abfs'' for usage\n');
        return
    end
end

% Count the number of files
nFiles = numel(fileNames);

%% Loop through all .abf files
abfParamsAllCell = cell(nFiles, 1);
dataAll = cell(nFiles, 1);
tVecAll = cell(nFiles, 1);
vVecsAll = cell(nFiles, 1);
iVecsAll = cell(nFiles, 1);
gVecsAll = cell(nFiles, 1);
dataReorderedAll = cell(nFiles, 1);
parfor iFile = 1:nFiles
%for iFile = 1:nFiles
    % Parse the abf file
    [abfParamsAllCell{iFile}, dataAll{iFile}, ...
        tVecAll{iFile}, vVecsAll{iFile}, ...
        iVecsAll{iFile}, gVecsAll{iFile}, dataReorderedAll{iFile}] = ...
        parse_abf(fileNames{iFile}, 'Verbose', verbose, ...
                    'UseOriginal', useOriginal, ...
                    'ExpMode', expMode, ...
                    'TimeUnits', timeUnits, ...
                    'ChannelTypes', channelTypes, ...
                    'ChannelUnits', channelUnits, ...
                    'ChannelLabels', channelLabels, ...
                    'IdentifyProtocols', identifyProtocols);
end

% Convert to a struct array
abfParamsAllStruct = [abfParamsAllCell{:}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

