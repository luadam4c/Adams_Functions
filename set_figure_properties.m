function fig = set_figure_properties (varargin)
%% Decides on the figure handle and sets figure properties
% Usage: fig = set_figure_properties (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       fig = set_figure_properties;
%       fig = set_figure_properties('Width', 200);
%       fig = set_figure_properties('Height', 300);
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
%                   - 'Position': figure position
%                   must be a 4-element positive integer vector
%                   default == get(0, 'defaultfigureposition')
%                   - 'Width': figure width
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (3)
%                   - 'Height': figure height
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (4)
%                   - Any other parameter-value pair for the figure() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_subplots.m
%       cd/plot_bar.m
%       cd/plot_frame.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-05-10 Created by Adam Lu
% 2019-08-23 Added 'FigExpansion' as an optional argument
% 2019-08-23 Renamed decide_on_fig_handle.m to set_figure_properties.m
% 2019-08-24 Now uses the default figure position
% 2019-09-04 Added 'Height', 'Width', 'Position' as optional arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figExpansionDefault = [];       % no figure expansion by default
positionDefault = [];           % set later
widthDefault = [];              % set later
heightDefault = [];             % set later

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
addParameter(iP, 'Position', positionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'Position must be a empty or a numeric vector!'));
addParameter(iP, 'Width', widthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'Height', heightDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figExpansion = iP.Results.FigExpansion;
position = iP.Results.Position;
width = iP.Results.Width;
height = iP.Results.Height;

% Keep unmatched arguments for the figure() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default figure position
if isempty(position)
    position = get(0, 'defaultfigureposition');

    % Modify figure position if requested
    if ~isempty(width)
        position(3) = width;
    end
    if ~isempty(height)
        position(4) = height;
    end
end

%% Do the job
% Decide on the figure handle
if ~isempty(figHandle)
    % Use the given figure
    fig = figure(figHandle, 'Position', position, otherArguments{:});
elseif ~isempty(figNumber)
    % Create and clear a new figure with given figure number
    fig = figure(figNumber, 'Position', position, otherArguments{:});

    % TODO: Make optional argument with different defaults
    clf(fig);
else
    % Get the current figure or create one if non-existent
    fig = gcf;

    % Set Figure object properties
    set(fig, 'Position', position, otherArguments{:});
end

% Expand figure position if requested
if ~isempty(figExpansion)
    fig = expand_figure_position(fig, figExpansion, position);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig = expand_figure_position (fig, expansionFactor, positionOld)
%% Expands or shrinks the figure position

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
