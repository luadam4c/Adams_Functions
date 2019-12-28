function [status, cmdOut] = compile_mod_files (varargin)
%% Compiles NEURON .mod files
% Usage: [status, cmdOut] = compile_mod_files (directory (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       status      - command status (0 is success, NaN is not run)
%                   specified as a nonnegative integer
%       cmdOut      - command output
%                   specified as a character vector
%
% Arguments:
%       directory   - (opt) directory containing .mod files
%                   must be a string scalar or a character vector
%                   default == pwd
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/all_files.m
%
% Used by:
%       cd/m3ha_network_launch.m
%       cd/m3ha_plot_figure03.m.m
%       cd/m3ha_rank_neurons.m
%       cd/run_neuron.m
%       /media/adamX/optimizer4compgabab/singleneuronfitting59.m

% File History:
% 2019-11-18 Adapted from m3ha_network_launch.m
% 


%% Hard-coded parameters
% TODO: Make optional argument
verbose = true;

%% Default values for optional arguments
directoryDefault = pwd;
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addOptional(iP, 'directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.directory;
% param1 = iP.Results.param1;

%% Preparation
% Initialize output
status = NaN;
cmdOut = '';

% Don't do anything if not running on a UNIX system
if ~isunix
    return
end

%% Do the job
% Display message
fprintf('Compiling NEURON .mod files in %s ...\n', directory);

% Test if there is a .mod file present in the home directory
[~, modPaths] = all_files('Directory', directory, 'Extension', 'mod');

% Change to the home directory and compile NEURON mod files
if isempty(modPaths)
    if verbose
        fprintf('Warning: There are no .mod files in %s!\n', directory);
    end
else
    cd(directory);
    [status, cmdOut] = unix('nrnivmodl');
end

% TODO: print status

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%