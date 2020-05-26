function [dydxVec, xVecNew] = compute_derivative_trace (yVec, xVec, varargin)
%% Computes the derivative trace dy/dx from x and y, using the midpoints of x as new x
% Usage: [dydxVec, xVecNew] = compute_derivative_trace (yVec, xVec, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       xVec = (0:100)';
%       yVec = sin(xVec);
%       [dydxVec, xVec1] = compute_derivative_trace(yVec, xVec);
%       [d2ydx2Vec, xVec2] = compute_derivative_trace(dydxVec, xVec1);
%
% Outputs:
%       dydxVec     - TODO: Description of dydxVec
%                   specified as a TODO
%       xVecNew     - TODO: Description of xVecNew
%                   specified as a TODO
%
% Arguments:
%       yVec        - y vector(s)
%                   must be a TODO
%       xVec        - x vector(s)
%                   must be a TODO
%       varargin    - 'FiltWidthSamples': moving average filter width in sample
%                   must be a numeric vector
%                   default == no filter
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/movingaveragefilter.m
%
% Used by:
%       cd/m3ha_find_decision_point.m
%       cd/m3ha_plot_simulated_traces.m

% File History:
% 2020-04-15 Created by Adam Lu
% 2020-05-19 Added 'FiltWidthSamples' as an optional argument

%% Hard-coded parameters

%% Default values for optional arguments
filtWidthSamplesDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'yVec');
addRequired(iP, 'xVec');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FiltWidthSamples', filtWidthSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, yVec, xVec, varargin{:});
filtWidthSamples = iP.Results.FiltWidthSamples;

%% Preparation
% Force as matrix
[xVec, yVec] = argfun(@force_matrix, xVec, yVec);

% Force as column vectors
[xVec, yVec] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonvectors', true), xVec, yVec);

%% Do the job
% Compute the new x points
xVecNew = (xVec(1:end-1, :) + xVec(2:end, :)) / 2;

% Compute the instantaneous derivative
dydxVec = diff(yVec) ./ diff(xVec);

%% Filter if requested
if ~isempty(filtWidthSamples)
    dydxVec = movingaveragefilter(dydxVec, filtWidthSamples);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%