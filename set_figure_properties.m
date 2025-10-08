function fig = set_figure_properties (varargin)
%% Decides on the figure handle and sets figure properties
% Usage: fig = set_figure_properties (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       fig = set_figure_properties;
%       fig = set_figure_properties('Width', 200);
%       fig = set_figure_properties('Units', 'centimeters', 'Width', 2, 'Height', 1);
%       fig = set_figure_properties('Height', 300);
%       fig = set_figure_properties('FigExpansion', [2, 2]);
%
% Outputs:
%       fig         - figure handle to use
%                   specified as a Figure object handle
%
% Arguments:
%       varargin    - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigExpansion': expansion factors for figure position
%                       Note: This occurs AFTER position is set
%                   must be a must be a positive scalar or 2-element vector
%                   default == []
%                   - 'ExpandFromDefault': whether to expand from figure 
%                                           position default
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true except when 'Position', 'Width' or 'Height'
%                               are set
%                   - 'Position': figure position
%                   must be a 4-element positive integer vector
%                   default == get(0, 'defaultfigureposition')
%                   - 'Width': figure width
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (3)
%                   - 'Height': figure height
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (4)
%                   - 'AdjustPosition': whether to adjust figure position 
%                                           so that it fits
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if 'FigExpansion', 'Width', 
%                               or 'Height' provided, but false otherwise
%                   - 'ClearFigure': whether to clear figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if 'FigNumber' provided 
%                               but false otherwise
%                   - 'AlwaysNew': whether to always create a new figure even if
%                                   figNumber is not passed in
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other properties for the Figure object
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compare_events_pre_post_stim.m
%       cd/create_subplots.m
%       cd/isemptycell.m
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/m3ha_fminsearch3.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_network_show_net.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_grouped_scatter.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_plot_violin.m
%       cd/m3ha_simulate_population.m
%       cd/parse_current_family.m
%       cd/parse_ipsc.m
%       cd/plot_autocorrelogram.m
%       cd/plot_bar.m
%       cd/plot_chevron.m
%       cd/plot_frame.m
%       cd/plot_histogram.m
%       cd/plot_measures.m
%       cd/plot_raster.m
%       cd/plot_relative_events.m
%       cd/plot_swd_histogram.m
%       cd/plot_traces.m
%       cd/plot_traces_spike2_mat.m
%       cd/plot_tuning_curve.m
%       cd/update_figure_for_corel.m
%       cd/virt_analyze_sniff_whisk.m
%       cd/virt_plot_log_decrement_jitter.m
%       cd/virt_plot_phase_response.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting.m
%       \Shared\Code\vIRt\virt_moore.m

% File History:
% 2019-05-10 Created by Adam Lu
% 2019-08-23 Added 'FigExpansion' as an optional argument
% 2019-08-23 Renamed decide_on_fig_handle.m to set_figure_properties.m
% 2019-08-24 Now uses the default figure position
% 2019-09-04 Added 'Height', 'Width', 'Position' as optional arguments
% 2019-09-06 Allowed 'FigExpansion' to be two elements
% 2019-09-06 Added 'AdjustPosition' and 'ClearFigure' as optional arguments
% 2019-09-08 Added 'AlwaysNew' as an optional argument
% 2019-09-12 Added 'ExpandFromDefault' as an optional argument
% 2019-11-17 Fixed 'AlwaysNew' when 'FigNumber' is provided
% TODO: Change axes outerPosition by default?

%% Hard-coded parameters

%% Default values for optional arguments
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figExpansionDefault = [];       % no figure expansion by default
expandFromDefaultDefault = [];  % set later
positionDefault = [];           % set later
widthDefault = [];              % set later
heightDefault = [];             % set later
adjustPositionDefault = [];     % set later
clearFigureDefault = [];        % set later
alwaysNewDefault = false;       % don't always create new figure by default

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
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'FigExpansion must be a empty or a numeric vector!'));
addParameter(iP, 'ExpandFromDefault', expandFromDefaultDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Position', positionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'Position must be a empty or a numeric vector!'));
addParameter(iP, 'Width', widthDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'Width must be a empty or a positive scalar!'));
addParameter(iP, 'Height', heightDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'Height must be a empty or a positive scalar!'));
addParameter(iP, 'AdjustPosition', adjustPositionDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ClearFigure', clearFigureDefault);
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figExpansion = iP.Results.FigExpansion;
expandFromDefault = iP.Results.ExpandFromDefault;
positionUser = iP.Results.Position;
width = iP.Results.Width;
height = iP.Results.Height;
adjustPosition = iP.Results.AdjustPosition;
clearFigureUser = iP.Results.ClearFigure;
alwaysNew = iP.Results.AlwaysNew;

% Keep unmatched name-value pairs for the Figure object
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% If alwaysNew is true, figHandle shouldn't be provided
if alwaysNew && ~isempty(figHandle)
    error('FigHandle can''t be provided if alwaysNew is true!');
end

% Set default expandFromDefault
if isempty(expandFromDefault)
    if ~isempty(figHandle) || ~isempty(positionUser) || ...
            ~isempty(width) || ~isempty(height)
        expandFromDefault = false;
    else
        expandFromDefault = true;
    end
end

% Decide whether to adjust figure position at the end
if isempty(adjustPosition)
    if ~isempty(width) || ~isempty(height) || ~isempty(figExpansion)
        adjustPosition = true;
    else
        adjustPosition = false;
    end
end

%% Do the job
% Decide on the figure handle
if ~isempty(figHandle)
    % Use the given figure
    fig = figure(figHandle);
elseif ~isempty(figNumber)
    % Close any figure opened with the same number
    if alwaysNew
        close(figure(figNumber));
    end

    % Create or go to the figure with given figure number
    fig = figure(figNumber);
elseif alwaysNew
    fig = figure;
else
    % Get the current figure or create one if non-existent
    fig = gcf;
end

% Decide whether to clear figure at the end
if isempty(clearFigureUser)
    if ~isempty(figNumber)
        clearFigure = true;
    else
        clearFigure = false;
    end
else
    clearFigure = clearFigureUser;
end

% Set other Figure object properties
if ~isemptycell(otherArguments)
    set(fig, otherArguments{:});
end

% Set figure position if requested
if ~isempty(positionUser)
    set(fig, 'Position', positionUser);
end

% Update figure width if requested
if ~isempty(width)
    positionNew = get(fig, 'Position');
    positionNew(3) = width;
    set(fig, 'Position', positionNew);
end

% Update figure height if requested
if ~isempty(height)
    positionNew = get(fig, 'Position');
    positionNew(4) = height;
    set(fig, 'Position', positionNew);
end

% Expand figure position if requested
if ~isempty(figExpansion)
    if expandFromDefault
        positionOld = get(0, 'defaultfigureposition');
    else
        positionOld = get(fig, 'Position');
    end

    expand_figure_position(fig, figExpansion, positionOld);
end

%% Adjust the figure position if needed
if adjustPosition
    adjust_figure_position(fig);
end

%% Clear the figure if requested
if clearFigure
    clf(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function expand_figure_position (fig, figExpansion, positionOld)
%% Expands or shrinks the figure position
% TODO: Pull out to its own function

% Force as a column vector
figExpansion = figExpansion(:);

% Make sure there are two elements
figExpansion = match_row_count(figExpansion, 2, ...
                                    'ExpansionMethod', 'repeat');

% Force as a row vector
figExpansion = transpose(figExpansion);

% Initialize as old position
if ~isempty(positionOld)
    positionNew = positionOld;
else
    positionNew = get(fig, 'Position');
end

% Compute a new figure length and width
positionNew(3:4) = positionOld(3:4) .* figExpansion;

% Compute the amount to shift starting points
positionShift = ((1 - figExpansion) ./ 2) .* positionOld(3:4);

% Compute a new figure starting points
positionNew(1:2) = positionOld(1:2) + positionShift;

% Set as new position
set(fig, 'Position', positionNew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figPositionNew = adjust_figure_position (fig)
%% Adjusts the figure position and return the new position
% TODO: Pull out to its own function

% Get the screen size
screenSize = get(0, 'ScreenSize');

% Get the current figure position
figPositionNow = get(fig, 'Position');

% Rename some positions
figLeft = figPositionNow(1);
figBottom = figPositionNow(2);
figHeight = figPositionNow(4);
figTop = figBottom + figHeight;

% Move the figure so that the top left is entirely on screen
if ~any(figPositionNow(3:4) > screenSize(3:4))
    % If the new figure size is not greater than the screen, use movegui()
    movegui(fig);

    % Get the new figure position
    figPositionNew = get(fig, 'Position');
else
    % Initialize the new figure position
    figPositionNew = figPositionNow;

    % If the top is out of screen, move to 
    if figLeft < screenSize(1) || figLeft > screenSize(1) + screenSize(3)
        figPositionNew(1) = screenSize(1);
    elseif figTop < screenSize(2) || figTop > screenSize(2) + screenSize(4)
        figPositionNew(2) = screenSize(2) + screenSize(4) - figHeight * 1.1;
    end

    % Update the new figure position
    set(fig, 'Position', figPositionNew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
OLD CODE:

if ~isempty(positionUser)
    positionOld = positionUser;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
