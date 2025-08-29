function lines = plot_error_bar (pValue, rLow, rHigh, varargin)
%% Plots error bar(s)
% Usage: lines = plot_error_bar (pValue, rLow, rHigh, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       lines       - Line objects plotted
%                       1 - connecting
%                       2 - upper limit
%                       3 - lower limit
%                   specified as a Line object handle array
% Arguments:
%       pValue      - parameter value of the error bar(s)
%                   must be a numeric vector
%       rLow        - lower readout value of the error bar(s)
%                   must be a numeric vector
%       rHigh       - upper readout value of the error bar(s)
%                   must be a numeric vector
%       varargin    - 'BarDirection': bar direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'vertical'   - vertical bars
%                       'horizontal' - horizontal bars
%                   default == 'vertical'
%                   - 'RelativeBarWidth': bar width relative to bar separation
%                   must be a nonnegative scalar
%                   default == 0.25
%                   - 'BarWidth': bar width(s)
%                   must be a nonnegative scalar
%                   default == 0.25 of bar separation
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       ~/Adams_Functions/argfun.m
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/force_column_vector.m
%       ~/Adams_Functions/plot_vertical_line.m
%       ~/Adams_Functions/plot_horizontal_line.m
%       ~/Adams_Functions/struct2arglist.m
%
% Used by:
%       ~/Adams_Functions/plot_bar.m
%       ~/Adams_Functions/plot_chevron.m

% File History:
% 2019-01-24 Created by Adam and Katerina
% 2019-05-10 Added 'BarDirection' as an optional argument
% 

%% Hard-coded parameters
validBarDirections = {'vertical', 'horizontal'};

%% Default values for optional arguments
barDirectionDefault = 'vertical';
relativeBarWidthDefault = 0.25;     % 0.25 of bar separation by default
barWidthDefault = [];               % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'pValue');
addRequired(iP, 'rLow');
addRequired(iP, 'rHigh');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BarDirection', barDirectionDefault, ...
    @(x) any(validatestring(x, validBarDirections)));
addParameter(iP, 'RelativeBarWidth', relativeBarWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'BarWidth', barWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, pValue, rLow, rHigh, varargin{:});
barDirection = validatestring(iP.Results.BarDirection, validBarDirections);
relativeBarWidth = iP.Results.RelativeBarWidth;
barWidth = iP.Results.BarWidth;

% Keep unmatched arguments for the line() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
if numel(pValue) ~= numel(rLow)
    disp("STOP NOW because number x values does not equal number of y Lows")
    lines = gobjects;
    return
elseif numel(pValue) ~= numel(rHigh)
    disp("STOP NOW because number x values does not equal number of y Highs")
    lines = gobjects;
    return
end

%% Preparation
% Compute bar separation
% TODO: Make this an optional parameter
barSeparation = min(diff(pValue));

% Set default bar width
if isempty(barWidth)
    % Set default
    barWidth = relativeBarWidth * barSeparation;
else
    % Don't overwrite default
end

% Count the number of lines
nLines = numel(pValue);

% Apply force_column_vector to each input argument
[pValue, rLow, rHigh] = argfun(@force_column_vector, pValue, rLow, rHigh);

% Compute the lower and upper parameter limits for the bars
pLow = pValue - 0.5 * barWidth;
pHigh = pValue + 0.5 * barWidth;

% Set readout limits
rLimits = transpose([rLow, rHigh]);

% Set parameter limits
pLimits = transpose([pLow, pHigh]);

% Preallocate an array for Line objects
lines = gobjects(nLines, 3);

%% Do the job
switch barDirection
    case 'vertical'
        % Plot the vertical line
        lines(:, 1) = ...
            plot_vertical_line(pValue, 'YLimits', rLimits, otherArguments{:});

        % Plot the upper limit horizontal line
        lines(:, 2) = ...
            plot_horizontal_line(rHigh, 'XLimits', pLimits, otherArguments{:});

        % Plot the lower limit horizontal line
        lines(:, 3) = ...
            plot_horizontal_line(rLow, 'XLimits', pLimits, otherArguments{:});
    case 'horizontal'
        % Plot the horizontal line
        lines(:, 1) = ...
            plot_horizontal_line(pValue, 'XLimits', rLimits, otherArguments{:});

        % Plot the upper limit vertical line
        lines(:, 2) = ...
            plot_vertical_line(rHigh, 'YLimits', pLimits, otherArguments{:});

        % Plot the lower limit vertical line
        lines(:, 3) = ...
            plot_vertical_line(rLow, 'YLimits', pLimits, otherArguments{:});
    otherwise
        error('barDirection unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

pValue = force_column_vector(pValue)
rLow = force_column_vector(rLow);
rHigh = force_column_vector(rHigh);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%