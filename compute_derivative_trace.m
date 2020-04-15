function [dydxVec, xVecNew] = compute_derivative_trace (xVec, yVec, varargin)
%% Computes the derivative trace dy/dx from x and y, using the midpoints of x as new x
% Usage: [dydxVec, xVecNew] = compute_derivative_trace (xVec, yVec, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       xVec = (0:100)';
%       [dydxVec, xVec1] = compute_derivative_trace(xVec, sin(xVec));
%       [d2ydx2Vec, xVec2] = compute_derivative_trace(xVec1, dydxVec);
%
% Outputs:
%       dydxVec     - TODO: Description of dydxVec
%                   specified as a TODO
%       xVecNew     - TODO: Description of xVecNew
%                   specified as a TODO
%
% Arguments:
%       xVec        - x vector(s)
%                   must be a TODO
%       yVec        - y vector(s)
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
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
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'xVec');
addRequired(iP, 'yVec');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, xVec, yVec, varargin{:});
% param1 = iP.Results.param1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%