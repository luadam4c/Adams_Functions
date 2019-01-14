function argList = struct2arglist (structure, varargin)
%% Converts a scalar structure to an argument list
% Usage: argList = struct2arglist (structure, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       argList     - argument list (parameter value pairs)
%                   specified as a cell array of character vectors
% Arguments:
%       structure   - structure with field names as parameters
%                   must be a scalar structure
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%
% Used by:
%       cd/annotation_in_plot.m
%       cd/compute_bins.m
%       cd/convert_to_rank.m
%       cd/medianfilter.m
%       cd/plot_grouped_histogram.m
%       cd/solve_function_at_value.m

% File History:
% 2018-12-28 Moved from annotation_in_plot.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'structure', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, structure, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Get all the parameter names
names = fieldnames(structure);

% Get all the parameter values
values = struct2cell(structure);

% Force as a column cell array
argList = force_column_cell(transpose([names, values]), 'ToLinearize', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%