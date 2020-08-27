function newValues = convert_units (oldValues, oldUnits, newUnits, varargin)
%% Converts numeric values from one units to another
% Usage: newValues = convert_units (oldValues, oldUnits, newUnits, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       convert_units(magic(3), 'mQ', 'nQ')
%       convert_units(magic(3), 'ms', 's')
%       convert_units(magic(3), 'hour', 'min')
%       convert_units(magic(3), 'Hz', 'rpm')
%
% Outputs:
%       newValues   - new numeric values
%                   specified as an array
%
% Arguments:
%       oldValues   - old numeric values
%                   must be an array
%       oldUnits    - old units
%                   must be a string scalar or a character vector
%       newUnits    - new Units
%                   must be a string scalar or a character vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/apply_to_all_cells.m
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_common_suffix.m
%
% Used by:
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_grouped_scatter.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_plot_violin.m
%       cd/minEASE.m
%       cd/parse_pleth_trace.m
%       cd/plot_raw_multiunit.m
%       cd/plot_relative_events.m
%       cd/plot_spectrogram_multiunit.m
%       cd/read_data_atf.m

% File History:
% 2020-01-22 Created by Adam Lu
% 2020-08-13 Now converts between all time units
% TODO: Convert between metric and non-metric units
% TODO: Use transform_vectors.m?

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'oldValues');
addRequired(iP, 'oldUnits', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'newUnits', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, oldValues, oldUnits, newUnits, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

% Check relationships between arguments
% TODO

%% Preparation
% Return old values if no computation necessary
if isempty(oldValues) || iscell(oldValues) && all(isemptycell(oldValues)) || ...
        strcmp(oldUnits, newUnits)
    newValues = oldValues;
    return
end

%% Do the job
% Extract the common suffix
commonSuffix = extract_common_suffix({oldUnits, newUnits}, 'Delimiter', '');

% Extract the distinct prefixes
if isempty(commonSuffix)
    oldPrefix = oldUnits;
    newPrefix = newUnits;
else
    oldPrefix = extractBefore(oldUnits, commonSuffix);
    newPrefix = extractBefore(newUnits, commonSuffix);
end

% Determine the old and new magnitudes
[oldMag, newMag] = argfun(@decide_on_magnitude, oldPrefix, newPrefix);

% Compute the appropriate scale factor
if ~isempty(oldMag) && ~isempty(newMag)
    scaleFactor = oldMag / newMag;
else
    error('Not implemented yet!');
end
    
% Return old values if no computation necessary
if scaleFactor == 1
    newValues = oldValues;
    return
end

% Multiply by the scaling factor
newValues = apply_to_all_cells(@(x) x .* scaleFactor, oldValues);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function magnitude = decide_on_magnitude (prefix)

%% Hard-coded constants
MIN_PER_S = 60;
HOUR_PER_S = 3600;
DAY_PER_S = 86400;
WEEK_PER_S = 604800;
YEAR_PER_S = 31536000;

%% Do the job
switch prefix
    case 'G'
        magnitude = 1e9;
    case 'M'
        magnitude = 1e6;
    case ''
        magnitude = 1;
    case 'd'
        magnitude = 1e-1;
    case 'c'
        magnitude = 1e-2;
    case 'm'
        magnitude = 1e-3;
    case 'u'
        magnitude = 1e-6;
    case 'n'
        magnitude = 1e-9;
    case 'p'
        magnitude = 1e-12;
    case 's'
        magnitude = 1;
    case 'min'
        magnitude = MIN_PER_S;
    case 'hour'
        magnitude = HOUR_PER_S;
    case 'day'
        magnitude = DAY_PER_S;
    case 'week'
        magnitude = WEEK_PER_S;
    case 'year'
        magnitude = YEAR_PER_S;
    case 'Hz'
        magnitude = 1;
    case 'rpm'
        magnitude = 1/60;
    otherwise
        error('prefix not implemented yet!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
