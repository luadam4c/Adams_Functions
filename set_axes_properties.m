function ax = set_axes_properties (varargin)
%% Decides on the axes handle and sets axes properties
% Usage: ax = set_axes_properties (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       ax = set_axes_properties;
%       ax = set_axes_properties('SubPlotNumber', [2, 3, 1]);
%
% Outputs:
%       ax         - axes handle to use
%                   specified as a Axes object handle
% Arguments:
%       varargin    - 'AxesHandle': axes handle for created axes
%                   must be a empty or a axes object handle
%                   default == []
%                   - 'SubPlotNumber': subplot numbers for created axes
%                   must be a 3-element positive integer vector
%                   default == []
%                   - 'Position': axes position
%                   must be a 4-element positive integer vector
%                   default == [0, 0, 1, 1] in normalized units
%                   - 'Width': axes width
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (3)
%                   - 'Height': axes height
%                   must be a positive scalar
%                   default == get(0, 'defaultfigureposition') (4)
%                   - Any other parameter-value pair for the axes() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isnumericvector.m
%       cd/ispositiveintegervector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_frame.m

% File History:
% 2019-09-04 Adaped from set_figure_properties.m


%% Hard-coded parameters

%% Default values for optional arguments
axHandleDefault = [];           % no existing axes by default
subPlotNumberDefault = [];      % no subplot number by default
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
addParameter(iP, 'AxesHandle', axHandleDefault);
addParameter(iP, 'SubPlotNumber', subPlotNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x), ...
                'SubPlotNumber must be a empty or a positive integer vector!'));
addParameter(iP, 'Position', positionDefault, ...
    @(x) assert(isempty(x) || isnumericvector(x), ...
                'Position must be a empty or a numeric vector!'));
addParameter(iP, 'Width', widthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'Height', heightDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
axHandle = iP.Results.AxesHandle;
subPlotNumber = iP.Results.SubPlotNumber;
positionUser = iP.Results.Position;
width = iP.Results.Width;
height = iP.Results.Height;

% Keep unmatched arguments for the axes() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set the default position to have no margins
% TODO: How to tell if a figure was newly created?
% if isempty(positionUser)
%     position = [0, 0, 1, 1];
% else
%     position = positionUser;
% end

%% Do the job
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

% Set other Axes object properties
set(ax, otherArguments{:});

% Modify Axes position if requested
% TODO: update_object_position.m
if ~isempty(positionUser)
    set(ax, 'Position', positionUser);
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
