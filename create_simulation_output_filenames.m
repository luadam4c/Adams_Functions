function outFilePaths = create_simulation_output_filenames (nSims, varargin)
%% Creates simulation output file names
% Usage: outFilePaths = create_simulation_output_filenames (nSims, varargin)
% Explanation:
%       TODO
% Example(s):
%       outFilePaths = create_simulation_output_filenames(nSims, ...
%                                   'OutFolder', outFolder, 'Prefix', prefix);
% Outputs:
%       outFilePaths - output file paths
%                   specified as a cell array of character vectors
% Arguments:    
%       nSims       - number of simulations
%                   must be a positive integer scalar
%       varargin    - 'Prefix': prefix to prepend to file names
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%
% Requires:
%       cd/create_labels_from_numbers.m
%       cd/force_string_end.m
%
% Used by:    
%       cd/m3ha_neuron_create_simulation_params.m

% File History:
% 2018-10-22 Created by Adam Lu
% 2018-12-17 Now uses create_labels_from_numbers.m
% TODO: Use force_string_start.m
% TODO (EASY): Make 'Suffix' an optional parameter with default ''
% TODO (EASY): Make 'Extension' an optional parameter with default '.out'
% 

%% Hard-coded parameters
suffix = '';
extension = '.out';

%% Default values for optional arguments
prefixDefault = '';             % prepend nothing to file names by default
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default

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
addRequired(iP, 'nSims', ...                  % number of simulations
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, nSims, varargin{:});
prefix = iP.Results.Prefix;
outFolder = iP.Results.OutFolder;

%% Preparation
% Make sure the prefix ends with a '_'
prefix = force_string_end(prefix, '_', 'OnlyIfNonempty', true);

% Make sure the suffix starts with a '_'
% prefix = force_string_start(suffix, '_', 'OnlyIfNonempty', true);

%% Do the job
% Generate a complete prefix including the path of the directory
prefixComplete = fullfile(outFolder, [prefix, 'sim']);

% Generate a complete suffix including the suffix and the file extension
suffixComplete = [suffix, extension];

% Create the file names iteratively
outFilePaths = create_labels_from_numbers(1:nSims, 'Prefix', prefixComplete, ...
                                        'Suffix', suffixComplete);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

outFilePaths = cell(nSims, 1);
parfor iSim = 1:nSims
    % Set output file path
    outFilePaths{iSim} = ...
        fullfile(outFolder, strcat(prefix, sprintf('sim%d.out', iSim)));
end

% Create an anonymous function that prints a file name from a number
myfun = @(x) fullfile(outFolder, strcat(prefix, sprintf('sim%d.out', x)));

% Label the file names from 1 to nSims
outFilePaths = arrayfun(myfun, transpose([1:nSims]), 'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
