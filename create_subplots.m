function [fig, ax] = create_subplots (varargin)
%% Creates subplots with maximal fit
% Usage: [fig, ax] = create_subplots (nPlotsOrNRows (opt), nColumns (opt), gridPositions (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [fig, ax] = create_subplots(4);
%       [fig, ax] = create_subplots(1, 4);
%       [fig, ax] = create_subplots(1, 1, 'FigNumber', 3);
%       [fig, ax] = create_subplots(1, 3, 'FigNumber', 4);
%       [fig, ax] = create_subplots(2, 3, 'FigNumber', 5);
%       [fig, ax] = create_subplots(2, 3, 2:3, 'FigNumber', 6);
%       [fig, ax] = create_subplots(2, 3, 'CenterPosition', [500, 500, 300, 200]);
%       [fig, ax] = create_subplots(6, 1, {1:3, 4, 5, 6}, 'CenterPosition', [500, 500, 400, 200]);
%
% Outputs:
%       fig         - figure handle
%                   specified as a figure handle
%       ax          - axes handle(s)
%                   specified as an array of axes handles
%
% Arguments:
%       nPlotsOrNRows   - (opt) number of plots 
%                           or number of rows if nVolumns provided
%                       must be a positive integer scalar
%                       default == 1
%       nColumns        - (opt) number of columns
%                       must be a positive integer scalar
%                       default == 1
%       gridPositions   - (opt) grid positions of subplots
%                       must be a positive integer vector 
%                           or a cell array of positive integer vectors
%                       default == num2cell(transpose(1:nGrids))
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
%                   - 'FigPosition': figure position
%                   must be a 4-element positive integer vector
%                   default == same as CenterPosition
%                   - 'FigWidth': figure width
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (3)
%                   - 'FigHeight': figure height
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (4)
%                   - 'ClearFigure': whether to clear figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'AlwaysNew': whether to always create a new figure even if
%                                   figNumber is not passed in
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'CenterPosition': position for the center subplot
%                   must be a 4-element positive integer vector
%                   default == get(0, 'defaultfigureposition')
%                   - 'AdjustPosition': whether to adjust figure position 
%                                           so that it fits
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for the subplot() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/isaninteger.m
%       cd/set_figure_properties.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/align_subplots.m
%       cd/crosscorr_profile.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_network_raster_plot.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_rank_neurons.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure08.m
%       cd/parse_current_family.m
%       cd/parse_multiunit.m
%       cd/parse_pleth_trace.m
%       cd/plot_table_parallel.m
%       cd/plot_relative_events.m
%       cd/plot_small_chevrons.m
%       cd/plot_spectrogram_multiunit.m
%       cd/plot_traces.m
%       cd/plot_traces_spike2_mat.m
%       cd/virt_golomb_generate_output.m

% File History:
% 2018-08-04 Created by Adam Lu
% 2019-09-06 Added 'FigPosition' and 'CenterPosition' as optional arguments
% 2019-09-06 Added gridPositions as an optional argument
% 2019-09-11 Added more figure properties as optional arguments
% 2020-02-06 Added 'ExpandFromDefault' as an optional argument
% 2020-04-19 Made the first argument nPlotsOrNRows 
% 2020-04-19 Made all arguments optional
% 2021-05-16 Added 'AdjustPosition' as an optional argument
% TODO: Added 'TransposeOrder' as an optional argument

%% Hard-coded parameters
horizontalDeadSpace = 0.25;     % relative dead space at the edges of figure 
                                %   but not in between subplots

%% Default values for optional arguments
nPlotOrNRowsDefault = 1;        % set later
nColumnsDefault = [];           % set later
gridPositionsDefault = [];      % set later
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figExpansionDefault = [];       % set later
expandFromDefaultDefault = [];  % set later
figPositionDefault = [];        % set later
figWidthDefault = [];           % set later
figHeightDefault = [];          % set later
clearFigureDefault = true;      % clear figure by default
alwaysNewDefault = false;       % don't always create new figure by default
centerPositionDefault = [];     % set later
adjustPositionDefault = true;   % adjust position to fit window by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'nPlotsOrNRows', nPlotOrNRowsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addOptional(iP, 'nColumns', nColumnsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addOptional(iP, 'gridPositions', gridPositionsDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x) || ...
                    iscellnumeric(x), ...
                ['gridPositions must be empty or a numeric vector ', ...
                    'or a cell array of positive integer vectors!']));


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
addParameter(iP, 'FigPosition', figPositionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'FigPosition must be a empty or a numeric vector!'));
addParameter(iP, 'FigWidth', figWidthDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'FigWidth must be a empty or a positive scalar!'));
addParameter(iP, 'FigHeight', figHeightDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'FigHeight must be a empty or a positive scalar!'));
addParameter(iP, 'ClearFigure', clearFigureDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CenterPosition', centerPositionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'Position must be a empty or a numeric vector!'));
addParameter(iP, 'AdjustPosition', adjustPositionDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
nPlotsOrNRows = iP.Results.nPlotsOrNRows;
nColumns = iP.Results.nColumns;
gridPositions = iP.Results.gridPositions;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figExpansion = iP.Results.FigExpansion;
expandFromDefault = iP.Results.ExpandFromDefault;
figPosition = iP.Results.FigPosition;
figWidth = iP.Results.FigWidth;
figHeight = iP.Results.FigHeight;
clearFigure = iP.Results.ClearFigure;
alwaysNew = iP.Results.AlwaysNew;
centerPosition = iP.Results.CenterPosition;
adjustPosition = iP.Results.AdjustPosition;

% Keep unmatched arguments for the subplot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on number of rows and columns
if isempty(nColumns)
    % The number of total numbers of subplots is given
    nPlots = nPlotsOrNRows;

    % Compute the square root
    sqrtNPlots = sqrt(nPlots);

    % Decide on the number of rows and columns
    if isaninteger(sqrtNPlots)
        nRows = sqrtNPlots;
        nColumns = sqrtNPlots;
    else
        nRows = nPlots;
        nColumns = 1;
    end
else
    % The number of rows is given
    nRows = nPlotsOrNRows;
end

% Count the number of grid positions
nGrids = nRows * nColumns;

% Decide on the center subplot position
if isempty(centerPosition)
    centerPosition = get(0, 'defaultfigureposition');
end

% Set default figure position and figure expansion factors
if isempty(figExpansion) && isempty(figPosition) && isempty(figHandle)
    % Start with the initial figure position
    figPosition = centerPosition;

    % Compute the horizontal expansion factor
    horizontalExpandFactor = (nColumns - (nColumns - 1) * horizontalDeadSpace);

    % Compute the vertical expansion factor
    verticalExpandFactor = nRows;

    % Compute the expansion factor
    figExpansion = [horizontalExpandFactor, verticalExpandFactor];
end

% Decide on the subplot gridPositions
if isempty(gridPositions)
    gridPositions = num2cell(transpose(1:nGrids));
else
    gridPositions = force_column_cell(gridPositions);
end

% Decide on the figure to plot on and set figure position
fig = set_figure_properties('FigHandle', figHandle, 'FigNumber', figNumber, ...
                    'FigExpansion', figExpansion, ...
                    'ExpandFromDefault', expandFromDefault, ...
                    'Position', figPosition, ...
                    'Width', figWidth, 'Height', figHeight, ...
                    'ClearFigure', clearFigure, 'AlwaysNew', alwaysNew, ...
                    'AdjustPosition', adjustPosition);

%% Create subplots
% Count the number of subplots
nSubPlots = numel(gridPositions);

% Find all axes in the figure
axOld = findall(fig, 'type', 'axes');

% If correct subplots already created, do nothing
if numel(axOld) == nSubPlots
    ax = axOld;
    return
end

% Create subplots
ax = arrayfun(@(x) create_one_subplot(x, gridPositions, ...
                                        nRows, nColumns, otherArguments), ...
                transpose(1:nSubPlots));

% If any subplot got deleted, recreate it
%   Note: For some reason, subplots sometime disappear 
while any(~isvalid(ax))
    ax = arrayfun(@(x) create_subplot_again(ax(x), x, gridPositions, ...
                                            nRows, nColumns, otherArguments), ...
                    transpose(1:nSubPlots));
end

% Decide whether future updates should not change the axes position
% TODO: Use update_figure_for_corel when needed instead?
% activePositionProperty = 'Position'
% set(ax, 'ActivePositionProperty', activePositionProperty);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function axThis = create_one_subplot (iSubPlot, gridPositions, ...
                                        nRows, nColumns, otherArguments)
%% Creates one subplot

% Get the grid positions for this subplot
gridPositionsThis = gridPositions{iSubPlot};

% Get the column numbers
iColumns = mod((gridPositionsThis - 1), nColumns) + 1;

% Get the row numbers
iRows = ceil(gridPositionsThis ./ nColumns);

% Get the minimum and maximum row and column numbers
minColumn = min(iColumns);
maxColumn = max(iColumns);
minRow = min(iRows);
maxRow = max(iRows);

% Compute the outer position for maximal fit
outerPositionThis = [(minColumn - 1)/nColumns, ...
                    (nRows - maxRow)/nRows, ...
                    (maxColumn - minColumn + 1)/nColumns, ...
                    (maxRow - minRow + 1)/nRows];

% Create subplot
axThis = subplot(nRows, nColumns, gridPositionsThis, ...
                otherArguments{:});

% Modify the outer position
set(axThis, 'OuterPosition', outerPositionThis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function axThis = create_subplot_again (axThis, iSubPlot, gridPositions, ...
                                        nRows, nColumns, otherArguments)
%% Creates a subplot again if it is not valid

if ~isvalid(axThis)
    % Recreate and save subplot in array
    axThis = create_one_subplot(iSubPlot, gridPositions, ...
                                    nRows, nColumns, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% 
figPosition = expand_figure_position(centerPosition, horizontalExpandFactor, verticalExpandFactor)

% Decide on the creating order
%   Note: For some reason, subplots sometime disappear 
if nRows > nColumns
    indSubplots = nSubPlots:-1:1;
else
    indSubplots = 1:nSubPlots;
end

% Initialize the axes array
ax = gobjects(nSubPlots, 1);
% Create subplots
for iSubPlot = 1:nSubPlots
    % Create and save subplot in array
    ax(iSubPlot) = 
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
