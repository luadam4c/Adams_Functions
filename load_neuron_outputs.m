function [outputs, fullPaths] = load_neuron_outputs (varargin)
%% Loads .out files created by NEURON into a cell array
% Usage: [outputs, fullPaths] = load_neuron_outputs (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       outputs     - outputs of NEURON simulations
%                   specified as a numeric array 
%                       or a cell array of numeric arrays
%
% Arguments:
%       varargin    - 'Directories': the name of the directory(ies) containing 
%                                   the .out files, e.g. '20161216'
%                   must be a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == pwd
%                   - 'FileNames': names of .out files to load
%                   must be empty, a character vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from pwd
%                   - 'Verbose': whether to output parsed results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveAfterLoad': whether to remove .out files 
%                                           after loading
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'tVecs': time vectors to match
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (none provided)
%                   - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TimeWindows' - time window(s) to restrict to
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%
% Requires:
%       cd/array_fun.m
%       cd/construct_and_check_fullpath.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/is_in_parallel.m
%       cd/match_format_vector_sets.m
%       cd/match_time_points.m
%
% Used by:    
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_simulate_population.m

% File History:
% 2018-10-23 Adapted from code in run_neuron_once_4compgabab.m
% 2018-10-31 Went back to using parfor for loading
% 2018-11-16 Fixed directories and allowed it to be a cell array TODO: fix all_files?
% 2020-01-01 Now uses array_fun.m
% 2020-01-31 Added 'ForceCellOutput' as an optional argument
% 2020-02-18 Added 'TimeWindows' as an optional argument
% 2020-04-24 Now loads and processes output within array_fun

%% Hard-coded parameters
outputExtension = '.out';

%% Default values for optional arguments
directoriesDefault = '';            % set later
fileNamesDefault = {};              % detect from pwd by default
verboseDefault = false;             % print to standard output by default
removeAfterLoadDefault = false;     % don't remove .out files by default
tVecsDefault = [];
forceCellOutputDefault = false; % don't force output as a cell array by default
timeWindowsDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directories', directoriesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveAfterLoad', removeAfterLoadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TimeWindows', timeWindowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['TimeWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, varargin{:});
directories = iP.Results.Directories;
fileNames = iP.Results.FileNames;
verbose = iP.Results.Verbose;
removeAfterLoad = iP.Results.RemoveAfterLoad;
tVecs = iP.Results.tVecs;
forceCellOutput = iP.Results.ForceCellOutput;
timeWindows = iP.Results.TimeWindows;

%% Preparation
% Decide on the files to use
if isempty(fileNames)
    % Find all .out files in the directories
    [~, fileNames] = all_files('Directory', directories, ...
                                'Extension', outputExtension, ...
                                'Verbose', verbose);

    % Return usage message if no .out files found
    if isempty(fileNames)
        fprintf('Type ''help %s'' for usage\n', mfilename);
        outputs = {};
        fullPaths = {};
        return
    end
end

% Force as a cell array
if ischar(fileNames)
    % Place in cell array
    fileNames = {fileNames};
elseif iscell(fileNames)
    % Always force output as a cell array in this case
    forceCellOutput = true;
end

% Construct full paths and check whether the files exist
%   TODO: Expand to accept optional Suffix', etc.
[fullPaths, pathExists] = ...
    construct_and_check_fullpath(fileNames, 'Directory', directories);

% Return if not all paths exist
if ~all(pathExists)
    fprintf('Some of the output paths do not exist!\n');
    outputs = {};
    fullPaths = {};
    return
end

% Match the number of time vectors and simulated outputs
[tVecs, fullPaths] = match_format_vector_sets(tVecs, fullPaths);

% Match the number of time windows and simulated outputs
[timeWindows, fullPaths] = match_format_vector_sets(timeWindows, fullPaths);

%% Load files
% Load the data saved by NEURON to a .out file into a cell array
outputs = array_fun(@(x, y, z) load_one_neuron_output(x, y, z), ...
                    fullPaths, tVecs, timeWindows, 'UniformOutput', false);

%% Outputs
% Don't output as cell if not necessary
if ~forceCellOutput && numel(outputs) == 1
    outputs = outputs{1};
end

%% Remove files
% Remove .out files created by NEURON if not to be saved
%   Note: Never use parfor here, so don't use array_fun 
if removeAfterLoad
    cellfun(@delete, fullPaths, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = load_one_neuron_output (fullPath, tVec, timeWindow)
%% Loads NEURON outputs from one file

% Load output
output = load(fullPath);

% If tVecs not empty, interpolate simulated data to match the time points
if ~isempty(tVec)
    % Interpolated simulated data
    output = match_time_points(output, tVec);
end

% Restrict output
if ~isempty(timeWindow)
    % Extract time vectors if needed
    if isempty(tVec)
        tVec = output(:, 1);
    end

    % Find window endpoints
    endPoints = find_window_endpoints(timeWindow, tVec);

    % Restrict to those endpoints
    output = extract_subvectors(output, 'EndPoints', endPoints);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of output files
nFiles = numel(fullPaths);

simDataNeuron = cell(nFiles, 1);
parfor iFile = 1:nFiles
    simDataNeuron{iFile} = load(fullPaths{iFile});
end

parfor iFile = 1:nFiles
    delete(fullPaths{iFile});
end

% 2018-11-01 The following is slower than parfor for large files
% Load the data saved by NEURON to a .out file into a cell array
outputs = cellfun(@load, fullPaths, 'UniformOutput', false);

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
