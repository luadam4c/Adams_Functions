function cellNames = m3ha_extract_cell_name (strs, varargin)
%% Extracts the cell name from strings but ignores anything before filesep
% Usage: cellNames = m3ha_extract_cell_name (strs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       m3ha_extract_cell_name({'D101310_13', 'D101310_15'})
%       m3ha_extract_cell_name({'ACB123456/D101310_13', 'ACB123456/D101310_15'})
%       m3ha_extract_cell_name({'D101310_13', 'D101310_15'}, 'ForceSingleOutput', true)
%
% Outputs:
%       cellNames    - cell names
%                   specified as a character vector
%
% Arguments:
%       strs        - strings
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'ForceSingleOutput': whether to force as a single output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%       cd/extract_substrings.m
%
% Used by:
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_simulate_ipsc_response.m

% File History:
% 2019-12-21 Created by Adam Lu
% 

%% Hard-coded parameters
cellNamePattern = '[A-Z][0-9]{6}';

%% Default values for optional arguments
forceSingleOutputDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'strs', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'forceSingleOutput', forceSingleOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, strs, varargin{:});
forceSingleOutput = iP.Results.forceSingleOutput;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Ignore anything before filesep
strs = extract_fileparts(strs, 'pathbase');
strs = extract_fileparts(strs, 'dirbase');

% Extract the cell names
cellNames = extract_substrings(strs, 'RegExp', cellNamePattern);

% Reduce to one cell name if requested
if forceSingleOutput
    % Extract the common cell name
    cellNames = extract_common_prefix(cellNames);
    
    % If doesn't exist, change it to 'mixed'
    if isempty(cellNames)
        cellNames = 'mixed';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%