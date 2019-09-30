function [edgesNew, isUpdated] = adjust_edges (edgesOld, varargin)
%% Update histogram bin edges according to specific parameters
% Usage: [edgesNew, isUpdated] = adjust_edges (edgesOld, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       adjust_edges(linspace(-5.3, 5.3, 12), 'FixedEdges', 0)
%       adjust_edges(linspace(-5.3, 5.3, 12), 'FixedEdges', [0, 3])
%       adjust_edges(linspace(-5.3, 5.3, 12), 'FixedEdges', [0, 3, -5])
%       adjust_edges(linspace(-5.3, 5.3, 11), 'FixedEdges', [0, 3, -5]) TODO
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
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/force_column_vector.m
%       cd/unique_custom.m
%
% Used by:
%       cd/compute_bins.m
%       cd/compute_grouped_histcounts.m

% File History:
% 2019-09-08 Created by Adam Lu
% 2019-09-08 Renamed as adjust_edges.m
% 2019-09-08 Now accepts up to three values for 'FixedEdges'
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
if ~isempty(fixedEdges) && ~all(ismember(fixedEdges, edgesOld))
    % Sort the fixed edges
    fixedEdges = unique_custom(fixedEdges, 'sorted', 'IgnoreNaN', true);

    % Force as column vectors
    [edgesOld, fixedEdges] = argfun(@force_column_vector, edgesOld, fixedEdges);

    % Add fixed edges to the edges to include
    edgesToInclude = sort([edgesOld; fixedEdges]);

    % Get the center fixed edge
    fixedCenter = extract_elements(fixedEdges, 'center');

    % Count the number of edges greater than the center fixed edge
    nEdgesRight = length(find(edgesToInclude > fixedCenter));

    % Count the number of edges less than the center fixed edge
    nEdgesLeft = length(find(edgesToInclude < fixedCenter));

    % Extract the average bin width from the old bins
    binWidthOld = nanmean(diff(edgesOld));

    % Compute the differences between fixed edges
    diffs = diff(fixedEdges);

    % Compute new bin widths if necessary
    if ~isempty(diffs) && ~all(mod(diffs, binWidthOld) == 0)
        % Get the extreme fixed edges
        fixedLeft = fixedEdges(1);
        fixedRight = fixedEdges(end);

        % Count the number of bins in between left and center
        nBinsLeftToCenter = round((fixedCenter - fixedLeft) / binWidthOld);

        % Count the number of bins in between left and center
        nBinsCenterToRight = round((fixedRight - fixedCenter) / binWidthOld);

        % Compute the bin width on either side
        if nBinsLeftToCenter ~= 0 && nBinsCenterToRight ~= 0
            binWidthLeft = (fixedCenter - fixedLeft) / nBinsLeftToCenter;
            binWidthRight = (fixedRight - fixedCenter) / nBinsCenterToRight;
        elseif nBinsLeftToCenter ~= 0 && nBinsCenterToRight == 0
            binWidthLeft = (fixedCenter - fixedLeft) / nBinsLeftToCenter;
            binWidthRight = binWidthLeft;
        elseif nBinsLeftToCenter == 0 && nBinsCenterToRight ~= 0
            binWidthRight = (fixedRight - fixedCenter) / nBinsCenterToRight;
            binWidthLeft = binWidthRight;
        else
            binWidthLeft = binWidthOld;
            binWidthRight = binWidthOld;
        end
    else
        binWidthLeft = binWidthOld;
        binWidthRight = binWidthOld;        
    end

    % Compute the new bin limits
    minEdge = fixedCenter - nEdgesLeft * binWidthLeft;
    maxEdge = fixedCenter + nEdgesRight * binWidthRight;

    % Compute the new bin edges
    edgesNew = [minEdge:binWidthLeft:fixedCenter, ...
                fixedCenter:binWidthRight:maxEdge];
    edgesNew = unique(edgesNew, 'sorted');

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

% Count the number of bins in between left and center
nBinsLeftToCenter = length(find(edgesToInclude >= fixedLeft & ...
                                edgesToInclude < fixedCenter));

% Count the number of bins in between left and center
nBinsCenterToRight = length(find(edgesToInclude > fixedCenter & ...
                                edgesToInclude <= fixedRight));

% Compute the bin width on either side
if nBinsLeftToCenter ~= 0
    binWidthLeft = (fixedCenter - fixedLeft) / nBinsLeftToCenter;
end
if nBinsCenterToRight ~= 0
    binWidthRight = (fixedRight - fixedCenter) / nBinsCenterToRight;
end
if nBinsLeftToCenter == 0
    binWidthLeft = binWidthRight;
end
if nBinsCenterToRight == 0
    binWidthRight = binWidthLeft;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
