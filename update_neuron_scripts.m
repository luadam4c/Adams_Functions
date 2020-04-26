function update_neuron_scripts (dirFrom, dirTo, varargin)
%% Updates NEURON scripts from one directory to another and change to the latter directory
% Usage: update_neuron_scripts (dirFrom, dirTo, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Arguments:
%       dirFrom     - directory containing scripts to update
%                   must be a string scalar or a character vector
%       dirTo       - directory to update
%                   must be a string scalar or a character vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for compile_mod_files()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/compile_mod_files.m
%
% Used by:
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_simulate_population.m

% File History:
% 2019-12-29 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'dirFrom', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'dirTo', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, dirFrom, dirTo, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the compile_mod_files() function
otherArguments = iP.Unmatched;

%% Do the job
% Make sure NEURON scripts are up to date
if isunix
    unix(sprintf('rsync -avhu %s/*.tem %s/*.mod %s/*.hoc %s', ...
                dirFrom, dirFrom, dirFrom, dirTo));
else
    % TODO
    error('Not implemented yet!');
end

% Compile .mod scripts and change to that directory
compile_mod_files(dirTo, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
