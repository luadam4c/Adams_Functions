function [edgesNew, isUpdated] = update_edges (edgesOld, varargin)
%% Update histogram bin edges according to specific parameters
% Usage: [edgesNew, isUpdated] = update_edges (edgesOld, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       edgesNew    - new edges
%                   specified as a column vector
%       isUpdated   - whether edges were updated
%                   specified as a logical scalar
%
% Arguments:
%       edgesOld    - original edges
%                   must be a numeric, logical, datetime or duration vector
%       varargin    - 'FixedEdges': numbers that must exist in bin edges
%                   must be a numeric, logical, datetime or duration vector
%                   default == []
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%
% Used by:
%       cd/compute_bins.m
%       cd/compute_grouped_histcounts.m

% File History:
% 2019-09-08 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
fixedEdgesDefault = [];

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
addRequired(iP, 'edgesOld', ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FixedEdges', fixedEdgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));


% Read from the Input Parser
parse(iP, edgesOld, varargin{:});
fixedEdges = iP.Results.FixedEdges;

%% Do the job
if ~all(ismember(fixedEdges, edgesOld))
    % Sort the fixed edges
    fixedEdges = sort(fixedEdges);

    % Get the center fixed edge
    centerEdge = extract_elements(fixedEdges, 'center');

    % Count the number of edges greater than the center fixed edge
    nEdgesRight = length(find(edgesOld > centerEdge));

    % Count the number of edges less than the center fixed edge
    nEdgesLeft = length(find(edgesOld < centerEdge));

    % Extract the average bin width
    binWidth = nanmean(diff(edgesOld));

    % Compute a new bin width if necessary
    diffs = diff(fixedEdges);
    if ~isempty(diffs) && ~all(mod(diffs, binWidth) == 0);
        % TODO
    end

    % Compute the new bin limits
    minEdge = centerEdge - nEdgesLeft * binWidth;
    maxEdge = centerEdge + nEdgesRight * binWidth;

    % Compute the new bin edges
    edgesNew = minEdge:binWidth:maxEdge;

    % Edges are updated
    isUpdated = true;
else
    % Compute the new bin edges
    edgesNew = edgesOld;

    % Edges are not updated
    isUpdated = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%