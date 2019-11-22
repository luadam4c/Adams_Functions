function [parsedParams, parsedData] = m3ha_parse_mat (matfiles, varargin)
%% Parses and loads a set of matfiles for the GAT blockade project
% Usage: [parsedParams, parsedData] = m3ha_parse_mat (matfiles, varargin)
%
% Outputs:
%       parsedParams    - a structure containing the following fields:
%                           nSamples
%                           siMs
%                           nSamples2
%                           siMs2
%       parsedData      - a structure containing the following fields:
%                           tvec0
%                           tvec2
%                           gvec0s
%                           ivec0s
%                           vvec0s
%                           gvec1s
%                           ivec1s
%                           vvec1s
%                           gvec2s
%                           ivec2s
%                           vvec2s
%                           vvec3s
%
% Arguments:
%       matfiles    - .mat file name(s)
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%       varargin    - 'LoadWindow': window (ms) to load
%                   must be within range of [timeBase, timeMax]
%                   default == [timeBase, timeMax]

%
% Requires:
%       cd/print_cellstr.m
%
% Used by:
%       cd/m3ha_dclampDataExtractor.m
%       cd/m3ha_dclampPassiveFitter.m

% File History:
% 2016-11-07 Moved from dclampDataExtractor.m
% 2018-10-03 Renamed load_matfiles_part -> m3ha_parse_mat
% 2018-10-03 Now places outputs in parsedParams, parsedData
% 2018-10-04 Now uses construct_and_check_fullpath.m
% 2018-10-15 Renamed ndps -> nSamples
% 2018-10-15 Added 'LoadWindow' as an optional argument

%% Default values for optional arguments
loadWindowDefault = [];            % don't restrict by default

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

% Add required inputs to an input Parser
addRequired(iP, 'matfiles', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['matfiles must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LoadWindow', loadWindowDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, matfiles, varargin{:});
loadWindow = iP.Results.LoadWindow;

%% Preparation
% Initialize outputs as empty arrays
parsedParams = [];
parsedData = [];

% Return if the argument is empty
if isempty(matfiles)
    fprintf('Nothing to parse!!\n\n');
    return
end

% Check whether matfiles exist
[fullPaths, fileExists] = ...
    construct_and_check_fullpath(matfiles, 'Extension', '.mat');

% If there are files missing, print them out and return
if ~all(fileExists)
    % Get all the full paths to the missing files
    missingFiles = fullPaths(~fileExists);

    % Print message and return
    fprintf('These files are missing:\n');
    print_cellstr(missingFiles, 'Delimiter', '\n', 'Prefix', '\t', ...
                                'OmitBraces', true, 'OmitQuotes', true);
    return
end

% Count the number of matfiles
if iscell(matfiles)
    nFiles = numel(matfiles);
else
    nFiles = 1;
end

% Force into a cell array
if ~iscell(matfiles)
    matfiles = {matfiles};
end

%% Load and parse
% Load first matfile
m = matfile(matfiles{1});

% Obtain time vectors from .mat data of first matfile
tvec0 = m.d_orig(:, 1);
tvec1 = m.d_mf(:, 1);
tvec2 = m.d_mfrs(:, 1);
tvec3 = m.d_mfmaf(:, 1);

% Compute the number of data points
nSamples = length(tvec0);
nSamples2 = length(tvec2);

% Compute the sampling intervals in milliseconds
siMs = tvec0(2) - tvec0(1);
siMs2 = tvec2(2) - tvec2(1);

% Get the time just before the start
timeBase = tvec0(1) - siMs;

% Get the maximum time
timeMax = tvec0(end);

% Decide on the indices to load
if isempty(loadWindow)
    indToLoad = 1:nSamples;
    indToLoad2 = 1:nSamples2;
elseif loadWindow(1) >= timeBase || loadWindow(2) <= timeMax
    % Find the start and end indices of each time vector that is within
    %   the load window
    idxStart = find(tvec0 >= loadWindow(1), 1);
    idxStart2 = find(tvec2 >= loadWindow(1), 1);
    idxEnd = find(tvec0 <= loadWindow(2), 1, 'last');
    idxEnd2 = find(tvec2 <= loadWindow(2), 1, 'last');

    % Construct indices to load
    indToLoad = idxStart:idxEnd;
    indToLoad2 = idxStart2:idxEnd2;    
else
    fprintf('Load window out of bounds!!\n');
    fprintf('   Must be within the range %s ~ %s ms\n', timeBase, timeMax);
    return
end

% Compute the number of data points to load
nSamplesToLoad = length(indToLoad);
nSamplesToLoad2 = length(indToLoad2);

% Restrict time vectors
tvec0 = tvec0(indToLoad);
tvec1 = tvec1(indToLoad);
tvec2 = tvec2(indToLoad2);
tvec3 = tvec3(indToLoad);

% Load data vectors for all matfiles
gvec0s = zeros(nSamplesToLoad, nFiles);
ivec0s = zeros(nSamplesToLoad, nFiles);
vvec0s = zeros(nSamplesToLoad, nFiles);
gvec1s = zeros(nSamplesToLoad, nFiles);
ivec1s = zeros(nSamplesToLoad, nFiles);
vvec1s = zeros(nSamplesToLoad, nFiles);
gvec2s = zeros(nSamplesToLoad2, nFiles);
ivec2s = zeros(nSamplesToLoad2, nFiles);
vvec2s = zeros(nSamplesToLoad2, nFiles);
gvec3s = zeros(nSamplesToLoad, nFiles);
ivec3s = zeros(nSamplesToLoad, nFiles);
vvec3s = zeros(nSamplesToLoad, nFiles);            
parfor iFile = 1:nFiles        % FOR each sweep
    % Get the path to this .mat file
    thisFile = matfiles{iFile};
    
    % Load this .mat file
    m = matfile(thisFile);

    % Extract vectors
    gvec0s(:, iFile) = m.d_orig(indToLoad, 2);
    ivec0s(:, iFile) = m.d_orig(indToLoad, 3);
    vvec0s(:, iFile) = m.d_orig(indToLoad, 4);
    gvec1s(:, iFile) = m.d_mf(indToLoad, 2);
    ivec1s(:, iFile) = m.d_mf(indToLoad, 3);
    vvec1s(:, iFile) = m.d_mf(indToLoad, 4);
    gvec2s(:, iFile) = m.d_mfrs(indToLoad2, 2);
    ivec2s(:, iFile) = m.d_mfrs(indToLoad2, 3);
    vvec2s(:, iFile) = m.d_mfrs(indToLoad2, 4);
    gvec3s(:, iFile) = m.d_mfmaf(indToLoad, 2);
    ivec3s(:, iFile) = m.d_mfmaf(indToLoad, 3);
    vvec3s(:, iFile) = m.d_mfmaf(indToLoad, 4);
end

%% Store outputs
% Store parameters
parsedParams.nSamples = nSamples;
parsedParams.siMs = siMs;
parsedParams.nSamples2 = nSamples2;
parsedParams.siMs2 = siMs2;
parsedParams.loadWindow = loadWindow;
parsedParams.nSamplesToLoad = nSamplesToLoad;
parsedParams.nSamplesToLoad2 = nSamplesToLoad2;

% Store vectors
parsedData.indToLoad = indToLoad;
parsedData.indToLoad2 = indToLoad2;
parsedData.tvec0 = tvec0;
parsedData.gvec0s = gvec0s;
parsedData.ivec0s = ivec0s;
parsedData.vvec0s = vvec0s;
parsedData.tvec1 = tvec1;
parsedData.gvec1s = gvec1s;
parsedData.ivec1s = ivec1s;
parsedData.vvec1s = vvec1s;
parsedData.tvec2 = tvec2;
parsedData.gvec2s = gvec2s;
parsedData.ivec2s = ivec2s;
parsedData.vvec2s = vvec2s;
parsedData.tvec3 = tvec3;
parsedData.gvec3s = gvec3s;
parsedData.ivec3s = ivec3s;
parsedData.vvec3s = vvec3s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[nSamples, siMs, tvec0, gvec0s, ivec0s, vvec0s, gvec1s, ivec1s, ...
    vvec1s, nSamples2, siMs2, tvec2, gvec2s, ivec2s, vvec2s, vvec3s] = ...
    m3ha_parse_mat (newinfolder, datfn, nswps)

filebase = [datfn, '_', num2str(swp)];
thisFile = fullfile(newinfolder, [filebase, '.mat']);
if ~exist(thisFile, 'file')
    error(['This mat file: ', thisFile, ' is missing!!']);
end

% Check whether matfiles exist
thisFile = fullfile(newinfolder, [datfn, '_1.mat']);
if exist(thisFile, 'file') ~= 2
    error(['This mat file: ', thisFile, ' is missing!!']);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%