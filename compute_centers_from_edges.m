function centers = compute_centers_from_edges (edges, varargin)
%% Computes bin centers from bin edges
% Usage: centers = compute_centers_from_edges (edges, varargin)
% Explanation:
%       TODO
% Example(s):
%       compute_centers_from_edges(1:5)
%       compute_centers_from_edges((1:5)')
% Outputs:
%       centers     - bin centers
%                   specified as a column vector
% Arguments:
%       edges       - bin edges
%                   must be accepted by the mean() function
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%
% Used by:
%       cd/compute_grouped_histcounts.m
%       cd/plot_grouped_histogram.m
%       cd/plot_histogram_with_outliers.m

% File History:
% 2019-01-15 Created by Adam Lu
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
addRequired(iP, 'edges');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, edges, varargin{:});
% param1 = iP.Results.param1;

%% Preparation
% Record whether edges is a row vector
if isrow(edges)
    isRow = true;
else
    isRow = false;
end

% Force as a column vector
edges = force_column_vector(edges);

%% Do the job
% Extract the left edges
lefEdges = edges(1:end-1);

% Extract the right edges
rightEdges = edges(2:end);

% Compute the mean of each pair of edges
centers = mean([lefEdges, rightEdges], 2);

%% Output
if isRow
    centers = force_row_vector(centers);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%