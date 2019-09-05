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
%                   - Any other properties for the Figure object
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
iP.KeepUnmatched = true;                        % allow extraneous options

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
positionUser = iP.Results.Position;
width = iP.Results.Width;
height = iP.Results.Height;

% Keep unmatched name-value pairs for the Figure object
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Decide on the figure handle
if ~isempty(figHandle)
    % Use the given figure
    fig = figure(figHandle);
elseif ~isempty(figNumber)
    % Create and clear a new figure with given figure number
    fig = figure(figNumber);
else
    % Get the current figure or create one if non-existent
    fig = gcf;
end

% Clear figure if requested
% TODO: Make optional argument with different defaults
if ~isempty(figNumber)
    clf(fig);
end

% Set other Figure object properties
set(fig, otherArguments{:});

% Modify figure position if requested
% TODO: update_object_position.m
if ~isempty(positionUser)
    set(fig, 'Position', positionUser);
end
if ~isempty(width)
    positionNew = get(fig, 'Position');
    positionNew(3) = width;
    set(fig, 'Position', positionNew);
end
if ~isempty(height)
    positionNew = get(fig, 'Position');
    positionNew(4) = height;
    set(fig, 'Position', positionNew);
end

% Expand figure position if requested
if ~isempty(figExpansion)
    if ~isempty(positionUser)
        positionOld = positionUser;
    elseif ~isempty(width) || ~isempty(height)
        positionOld = get(fig, 'Position');
    else
        positionOld = get(0, 'defaultfigureposition');
    end

    fig = expand_figure_position(fig, figExpansion, positionOld);
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

% TODO: Make the following adjust_figure_position.m
%%
% Get the screen size
screenSize = get(0, 'ScreenSize');

% If the new figure size is not greater than the screen, use movegui()
if ~any(positionNew(3:4) > screenSize(3:4))
    % Move the figure to entirely on screen
    movegui(fig)
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
