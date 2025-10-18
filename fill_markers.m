function fill_markers (varargin)
%% Fills markers if any for the current axes
% Usage: fill_markers (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       fill_markers;
%
% Arguments:
%       varargin    - 'AxesHandle': axes handle for created axes
%                   must be a empty or a axes object handle
%                   default == []
%                   
% Requires:
%       cd/set_axes_properties.m
%
% Used by:
%       cd/plot_tuning_curve.m

% File History:
% 2019-10-03 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
axHandleDefault = [];           % no existing axes by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, varargin{:});
axHandle = iP.Results.AxesHandle;

%% Preparation
% Decide on the axes
ax = set_axes_properties('AxesHandle', axHandle);

%% Do the job
% Find all Line objects
lines = findobj(ax, 'Type', 'Line');

% Fill markers with the same color as the marker edge
arrayfun(@fill_markers_for_one_line, lines);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fill_markers_for_one_line(lineObject)
%% Fill markers for one line object

% Get the marker edge color
markerEdgeColor = get(lineObject, 'MarkerEdgeColor');

% If set to 'auto', get the line color
if ischar(markerEdgeColor) && strcmpi(markerEdgeColor, 'auto')
    markerEdgeColor = get(lineObject, 'Color');
end

% Set the marker face color to be the same as the edge
set(lineObject, 'MarkerFaceColor', markerEdgeColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%