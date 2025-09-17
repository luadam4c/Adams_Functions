function ax = set_axes_properties (varargin)
%% Decides on the axes handle and sets axes properties
% Usage: ax = set_axes_properties (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       ax = set_axes_properties;
%       ax = set_axes_properties('AxesHandle', gca);
%       ax = set_axes_properties('SubPlotNumber', [2, 3, 1]);
%       ax = set_axes_properties('AxesCoverage', 90);
%
% Outputs:
%       ax         - axes handle to use
%                   specified as a Axes object handle
%
% Arguments:
%       varargin    - 'AxesHandle': axes handle for created axes
%                   must be a empty or a axes object handle
%                   default == []
%                   - 'SubPlotNumber': subplot numbers for created axes
%                   must be a 3-element positive integer vector
%                   default == []
%                   - 'OuterPosition': axes outer position
%                   must be a 4-element positive integer vector
%                   default == [] (not set)
%                   - 'Position': axes position
%                   must be a 4-element positive integer vector
%                   default == [] (not set)
%                   - 'Width': axes width
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (3)
%                   - 'Height': axes height
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (4)
%                   - 'AxesCoverage': percent of axes position 
%                                       relative to outerPosition
%                   must be a nonnegative scalar or 2-element vector
%                   default == [] (not changed from default)
%                   - 'FigHandle': figure handle to look for axes
%                   must be a empty or a figure object handle
%                   default == []
%                   - Any other parameter-value pair for the axes() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isemptycell.m
%       cd/isnumericvector.m
%       cd/ispositiveintegervector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/fill_markers.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/plot_chevron.m
%       cd/plot_correlation_coefficient.m
%       cd/plot_frame.m
%       cd/plot_grouped_jitter.m
%       cd/plot_grouped_scatter.m
%       cd/plot_raster.m
%       cd/plot_traces.m
%       cd/plot_test_result.m
%       cd/plot_tuning_curve.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m

% File History:
% 2019-09-04 Adaped from set_figure_properties.m
% 2019-09-19 Now defaults TickDir to 'out'
% 2019-11-05 Added 'FigHandle' as an optional argument

%% Hard-coded parameters

%% Default values for optional arguments
axHandleDefault = [];           % gca by default
subPlotNumberDefault = [];      % no subplot number by default
outerPositionDefault = [];      % set later
positionDefault = [];           % set later
widthDefault = [];              % set later
heightDefault = [];             % set later
axesCoverageDefault = [];       % set later
figHandleDefault = [];          % no existing figure by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AxesHandle', axHandleDefault);
addParameter(iP, 'SubPlotNumber', subPlotNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x), ...
                'SubPlotNumber must be a empty or a positive integer vector!'));
addParameter(iP, 'OuterPosition', outerPositionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'OuterPosition must be a empty or a numeric vector!'));
addParameter(iP, 'Position', positionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'Position must be a empty or a numeric vector!'));
addParameter(iP, 'Width', widthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'Height', heightDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'AxesCoverage', axesCoverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'vector'}));
addParameter(iP, 'FigHandle', figHandleDefault);

% Read from the Input Parser
parse(iP, varargin{:});
axHandle = iP.Results.AxesHandle;
subPlotNumber = iP.Results.SubPlotNumber;
outerPosition = iP.Results.OuterPosition;
position = iP.Results.Position;
width = iP.Results.Width;
height = iP.Results.Height;
axesCoverage = iP.Results.AxesCoverage;
figHandle = iP.Results.FigHandle;

% Keep unmatched arguments for the axes() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation

%% Do the job
% Go to the designated figure
if ~isempty(figHandle)
    figure(figHandle);
end

% Decide on the axes handle
if ~isempty(axHandle)
    % Use the given axes
    axes(axHandle);

    % Return the handle
    ax = gca;
elseif ~isempty(subPlotNumber)
    % Put numbers in a cell array
    subPlotNumberCell = num2cell(subPlotNumber);

    % Create and clear a new axes with given subplot number
    ax = subplot(subPlotNumberCell{:});
else
    % Get the current axes or create one if non-existent
    ax = gca;
end

% Plot with tick direction outward
% TODO: Make optional argument with this default
set(ax, 'TickDir', 'out');
set(ax, 'TickDirMode', 'manual');

% Set other Axes object properties
if ~isemptycell(otherArguments)
    set(ax, otherArguments{:});
end

% Modify Axes position if requested
% TODO: update_object_position.m
if ~isempty(outerPosition)
    set(ax, 'OuterPosition', outerPosition);
end

% If no margins, make the position the same as outer position
if ~isempty(axesCoverage)
    outerPosition = get(ax, 'OuterPosition');
    position(1:2) = outerPosition(1:2) + ...
                        outerPosition(3:4) .* (100 - axesCoverage)/200;
    position(3:4) = outerPosition(3:4) .* axesCoverage/100;
end
if ~isempty(position)
    set(ax, 'Position', position);
end
if ~isempty(width)
    positionNew = get(ax, 'Position');
    positionNew(3) = width;
    set(ax, 'Position', positionNew);
end
if ~isempty(height)
    positionNew = get(ax, 'Position');
    positionNew(4) = height;
    set(ax, 'Position', positionNew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
