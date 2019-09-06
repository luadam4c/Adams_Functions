function [fig, ax] = create_subplots (nRows, nColumns, varargin)
%% Creates subplots with maximal fit
% Usage: [fig, ax] = create_subplots (nRows, nColumns, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [fig, ax] = create_subplots(1, 1, 'FigNumber', 3);
%       [fig, ax] = create_subplots(1, 3, 'FigNumber', 4);
%       [fig, ax] = create_subplots(2, 3, 'FigNumber', 5);
%       [fig, ax] = create_subplots(2, 3, 2:3, 'FigNumber', 6);
%       [fig, ax] = create_subplots(2, 3, 'CenterPosition', [500, 500, 300, 200]);
%
% Outputs:
%       fig         - figure handle
%                   specified as a figure handle
%       ax          - axes handle(s)
%                   specified as an array of axes handles
% Arguments:
%       nRows       - number of rows
%                   must be a positive integer scalar
%       nColumns    - number of columns
%                   must be a positive integer scalar
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
%                   - 'FigPosition': figure position
%                   must be a 4-element positive integer vector
%                   default == expanded from CenterPosition
%                   - 'CenterPosition': position for the center subplot
%                   must be a 4-element positive integer vector
%                   default == get(0, 'defaultfigureposition')
%                   - Any other parameter-value pair for the subplot() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/set_figure_properties.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/parse_multiunit.m

% File History:
% 2018-08-04 Created by Adam Lu
% 2018-09-06 Added 'FigPosition' and 'CenterPosition' as optional arguments
% 2018-09-06 Added gridPositions as an optional argument
% 

%% Hard-coded parameters
horizontalDeadSpace = 0.25;     % relative dead space at the edges of figure 
                                %   but not in between subplots

%% Default values for optional arguments
gridPositionsDefault = [];      % set later
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figPositionDefault = [];        % set later
centerPositionDefault = [];     % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'nRows', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addRequired(iP, 'nColumns', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Add optional inputs to the Input Parser
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
addParameter(iP, 'FigPosition', figPositionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'Position must be a empty or a numeric vector!'));
addParameter(iP, 'CenterPosition', centerPositionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'Position must be a empty or a numeric vector!'));

% Read from the Input Parser
parse(iP, nRows, nColumns, varargin{:});
gridPositions = iP.Results.gridPositions;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figPosition = iP.Results.FigPosition;
centerPosition = iP.Results.CenterPosition;

% Keep unmatched arguments for the subplot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Count the number of grid positions
nGrids = nRows * nColumns;

% Decide on the center subplot position
if isempty(centerPosition)
    centerPosition = get(0, 'defaultfigureposition');
end

% Decide on the figure position
if isempty(figPosition)
    % TODO: expand_figure_position.m
    % Compute the horizontal expansion factor
    horizontalExpandFactor = (nColumns - (nColumns - 1) * horizontalDeadSpace);

    % Compute the vertical expansion factor
    verticalExpandFactor = nRows;

    % Copy original position
    figPosition = centerPosition;

    % Modify new position
    figPosition(1) = centerPosition(1) - centerPosition(3);
    figPosition(2) = centerPosition(2) - centerPosition(4);
    figPosition(3) = horizontalExpandFactor * centerPosition(3);
    figPosition(4) = verticalExpandFactor * centerPosition(4);
end

% Decide on the subplot gridPositions
if isempty(gridPositions)
    gridPositions = num2cell(transpose(1:nGrids));
else
    gridPositions = force_column_cell(gridPositions);
end

% Decide on the figure to plot on and set figure position
fig = set_figure_properties('FigHandle', figHandle, 'FigNumber', figNumber, ...
                            'Position', figPosition);

% TODO: Make the following adjust_figure_position.m
%%
% Get the screen size
screenSize = get(0, 'ScreenSize');

% Get the current figure position
figPositionNow = get(fig, 'Position');

% If the new figure size is not greater than the screen, use movegui()
if ~any(figPositionNow(3:4) > screenSize(3:4))
    % Move the figure to entirely on screen
    movegui(fig);
end
%%

%% Create subplots
% Count the number of subplots
nSubPlots = numel(gridPositions);

% Initialize the axes array
ax = gobjects(nSubPlots, 1);

% Save each subplot
for iSubPlot = 1:nSubPlots
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

    % Create subplot
    axThis = subplot(nRows, nColumns, gridPositionsThis, ...
                    otherArguments{:});

    % Modify the outer position for maximal fit
    axThis.OuterPosition = [(minColumn - 1)/nColumns, ...
                            (nRows - maxRow)/nRows, ...
                            (maxColumn - minColumn + 1)/nColumns, ...
                            (maxRow - minRow + 1)/nRows];

    % Save subplot in array
    ax(iSubPlot) = axThis;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
