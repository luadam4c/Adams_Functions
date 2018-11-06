function varargout = argfun (myFunction, varargin)
%% Applies a function to each input argument
% Usage: varargout = argfun (myFunction, varargin)
% Explanation:
%       This function applies the first argument (a function)
%           to each of the rest of the arguments
%       The first argument becomes the first output,
%           the second argument becomes the second output, and so on ...
% Example(s):
%       [a, b] = argfun(@sum, 1:10, magic(3))
% Outputs:
%       varargout   - outputs of the function applied to each input argument
% Arguments:    
%       myFunction  - a custom function
%                   must be a function handle
%       varargin    - input arguments
%
% Used by:
%       cd/compute_rms_error.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/construct_fullpath.m
%       cd/find_window_endpoints.m
%       cd/m3ha_create_initial_neuronparams.m
%       cd/m3ha_plot_individual_traces.m
%       cd/m3ha_run_neuron_once.m
%       cd/match_format_vectors.m
%       cd/match_reciprocals.m
%       cd/plot_cfit_pulse_response.m

% File History:
% 2019-10-25 Created by Adam Lu
% 

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
addRequired(iP, 'myFunction', ...                  % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));

% Read from the Input Parser
parse(iP, myFunction);

%% Preparation
% Count the number of inputs (number of arguments excluding the function handle)
nInputs = nargin - 1;

% If there are more outputs requested than inputs, return error
if nargout > nInputs
    error(['Cannot request more outputs than provided inputs, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
% Generate field names for the input arguments
myFieldNames = arrayfun(@(x) ['Arg', num2str(x)], 1:nInputs, ...
                        'UniformOutput', false);

% Place all arguments in an input structure
%   Note: varargin is a row cell array
myStructInputs = cell2struct(varargin, myFieldNames, 2);

% Pass the function to all fields in the input structure
myStructOutputs = structfun(myFunction, myStructInputs, ...
                            'UniformOutput', false);

% Extract all fields from the structure
varargout = struct2cell(myStructOutputs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%