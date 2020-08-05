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
%       varargin    - 'FromBaseName': whether to restrict to the base name
%                                       (remove everything before last filesep
%                                        and after last .)
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for extract_substrings()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_substrings.m
%
% Used by:
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_simulate_population.m

% File History:
% 2019-12-21 Created by Adam Lu
% 2020-05-06 Added cellNamePatternAlt

%% Hard-coded parameters
cellNamePattern = '[A-Z][0-9]{6}';
cellNamePatternAlt = 'hetero[0-9]*seed[0-9]*';

%% Default values for optional arguments
fromBaseNameDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'strs', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FromBaseName', fromBaseNameDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, strs, varargin{:});
fromBaseName = iP.Results.FromBaseName;

% Keep unmatched arguments for the extract_substrings() function
otherArguments = iP.Unmatched;

%% Do the job
% Extract the cell names
cellNames = extract_substrings(strs, 'FromBaseName', fromBaseName, ...
                                'RegExp', cellNamePattern, otherArguments);

% Determine if there is any empty string
toReplace = isemptycell(cellNames);

% Replace with alternative cell name string
if any(toReplace)
    cellNamesAlt = extract_substrings(strs, 'FromBaseName', fromBaseName, ...
                                    'RegExp', cellNamePatternAlt, otherArguments);
    cellNames(toReplace) = cellNamesAlt(toReplace);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
