function fig = set_figure_properties (varargin)
%% Decides on the figure handle and sets figure properties
% Usage: fig = set_figure_properties (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       fig = set_figure_properties;
%       fig = set_figure_properties('FigExpansion', 2);
%
% Outputs:
%       fig         - figure handle to use
%                   specified as a Figure object handle
% Arguments:
%       varargin    - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigExpansion': expansion factor for figure position
%                   must be a positive scalar
%                   default == []
%                   - Any other parameter-value pair for the figure() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_subplots.m
%       cd/plot_bar.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-05-10 Created by Adam Lu
% 2019-08-23 Added 'FigExpansion' as an optional argument
% 2019-08-23 Renamed to set_figure_properties.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figExpansionDefault = [];       % no figure expansion by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'FigExpansion', figExpansionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figExpansion = iP.Results.FigExpansion;

% Keep unmatched arguments for the figure() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Decide on the figure handle
if ~isempty(figHandle)
    % Use the given figure
    fig = figure(figHandle, otherArguments{:});
elseif ~isempty(figNumber)
    % Create and clear a new figure with given figure number
    fig = figure(figNumber, otherArguments{:});
    clf(fig);
else
    % Get the current figure or create one if non-existent
    fig = gcf;
end

% Expand the figure position if requested
if ~isempty(figExpansion)
    fig = expand_figure_position(fig, figExpansion);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig = expand_figure_position (fig, expansionFactor)
%% Expands or shrinks the figure position

% Get the old figure position
positionOld = get(fig, 'Position');

% Initialize as old position
positionNew = positionOld;

% Compute a new figure length and width
positionNew(3:4) = positionOld(3:4) * expansionFactor;

% Compute the amount to shift starting points
positionShift = ((1 - expansionFactor) / 2) * positionOld(3:4);

% Compute a new figure starting points
positionNew(1:2) = positionOld(1:2) + positionShift;

% Set as new position
set(fig, 'Position', positionNew);

% Get the screen size
screenSize = get(0, 'ScreenSize');

% If the new figure size is not greater than the screen, use movegui()
if ~any(positionNew(3:4) > screenSize(3:4))
    % Move the figure to entirely on screen
    movegui(fig)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
