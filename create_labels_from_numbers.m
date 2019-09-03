function labels = create_labels_from_numbers (numbers, varargin)
%% Creates a cell array of labels from an array of numbers with an optional prefix or suffix
% Usage: labels = create_labels_from_numbers (numbers, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       labels = create_labels_from_numbers([3, 7, 5])
%       labels = create_labels_from_numbers(1:3, 'Suffix', ' Mississippi')
%       labels = create_labels_from_numbers(1:3, 'Prefix', 'Husband #')
%       labels = create_labels_from_numbers(1:3, 'Prefix', 'Katie ')
%       labels = create_labels_from_numbers(1:3, 'Prefix', 'Make ', 'Suffix', ' Wish')
%
% Outputs:
%       labels     - labels created
%                   specified as a cell array of character vectors
% Arguments:
%       numbers     - numbers
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'ForceColumnOutput': whether to force output as a column
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Prefix': string to place before each number
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'Suffix': string to place after each number
%                   must be a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for the convert_to_char() function
%
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/convert_to_char.m
%
% Used by:
%       cd/atfwrite.m
%       cd/combine_variables_across_tables.m
%       cd/compute_all_pulse_responses.m
%       cd/create_default_grouping.m
%       cd/create_row_labels.m
%       cd/create_simulation_output_filenames.m
%       cd/m3ha_xolotl_plot.m
%       cd/plot_histogram.m
%       cd/parse_lts.m
%       cd/plot_raster.m
%       cd/plot_struct.m
%       cd/plot_swd_raster.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%       cd/renamevars.m
%       cd/save_all_zooms.m

% File History:
% 2018-12-17 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
forceColumnOutputDefault = true;    % force output as a column by default
prefixDefault = '';     % no string to place before each number by default
suffixDefault = '';     % no string to place after each number by default

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
addRequired(iP, 'numbers', ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceColumnOutput', forceColumnOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, numbers, varargin{:});
forceColumnOutput = iP.Results.ForceColumnOutput;
prefix = iP.Results.Prefix;
suffix = iP.Results.Suffix;

% Keep unmatched arguments for the convert_to_char() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force the numbers as a column vector if requested
if forceColumnOutput
    numbers = numbers(:);
end

%% Do the job
% Create the labels
labels = arrayfun(@(x) [prefix, convert_to_char(x, otherArguments{:}), suffix], ...
                    numbers, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

addRequired(iP, 'nLabels', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

labels = arrayfun(@(x) [prefix, num2str(x), suffix], transpose(1:nLabels), ...
                        'UniformOutput', false);

% strcat omits spaces in some cases
labels = arrayfun(@(x) strcat(prefix, convert_to_char(x), suffix), ...
                    numbers, 'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
