function indZeros = find_zeros (vecs, varargin)
%% Find the indices where the vectors are closest to zero
% Usage: indZeros = find_zeros (vecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       indZeros    - TODO: Description of indZeros
%                   specified as a TODO
%
% Arguments:
%       vecs        -  TODO: Description of vecs
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       /TODO:dir/TODO:file
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m

% File History:
% 2020-04-15 Created by Adam Lu
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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vecs');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vecs, varargin{:});
% param1 = iP.Results.param1;

%% Preparation
% Force as column vectors
vecs = force_column_vector(vecs);

%% Do the job
if iscell(vecs)
    indZeros = cellfun(@find_zeros_helper, vecs, 'UniformOutput', false);
else
    indZeros = find_zeros_helper(vecs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indZeros = find_zeros_helper(vec)

% Find the left and right of each consecutive pair of sample points
vecLeft = vec(1:end-1);
vecRight = vec(2:end);

% Find the indices with vec values closest to zero
indLeftZeros = find(vecLeft .* vecRight <= 0);
indRightZeros = indLeftZeros + 1;

% Choose the index closest to zero
indZeros = arrayfun(@(a, b) choose_closest_to_value(vec, a, b, 0), ...
                    indLeftZeros, indRightZeros);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = choose_closest_to_value (vec, idx1, idx2, value)

if abs(vec(idx1) - value) <= abs(vec(idx2) - value)
    idx = idx1;
else
    idx = idx2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%