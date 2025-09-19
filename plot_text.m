function [htext] = plot_text (textStr, varargin)
%% Plots a text string on an axes
% Usage: [htext] = plot_text (textStr, varargin)
% Explanation:
%       This function plots a text string on a specified or current axes.
%       If a cell array of strings is provided, it concatenates them with a
%       delimiter. The position is determined by a location string and a margin.
%
% Example(s):
%       figure;
%       plot(1:10, 1:10);
%       plot_text('Sample Text');
%       plot_text({'Multiple', 'Lines'}, 'TextLocation', 'bottomright', 'Color', 'b');
%       plot_text('Custom Alignment', 'TextLocation', 'topright', 'HorizontalAlignment', 'center');
%       plot_text('Data Units', 'Units', 'data', 'TextLocation', 'bottomleft', 'TextMargin', 1); % Margin is now in data units
%
% Outputs:
%       htext       - Handle to the text object created
%                   specified as a Text object
%
% Arguments:
%       textStr     - The text string(s) to plot
%                   must be a character vector, a string scalar,
%                   or a cell array of character vectors/strings
%       varargin    - 'Delimiter': Delimiter for joining cell arrays of strings
%                   must be a character vector or string scalar
%                   default == '\n' (newline)
%                   - 'Units': coordinate system for text positioning
%                   must be an unambiguous, case-insensitive match to one of:
%                       'normalized'  - relative to the axes (default)
%                       'data'        - uses the data coordinates of the plot
%                       'pixels'      - uses pixel units
%                       'inches'      - uses inches
%                       'centimeters' - uses centimeters
%                       'points'      - uses points
%                   default == 'normalized'
%                   - 'TextLocation': location for the text
%                   must be an unambiguous, case-insensitive match to one of:
%                       'topleft'     - places text in the top-left corner
%                       'topright'    - places text in the top-right corner
%                       'bottomleft'  - places text in the bottom-left corner
%                       'bottomright' - places text in the bottom-right corner
%                   default == 'topleft'
%                   - 'TextMargin': Margin for the text location
%                   must be a numeric scalar
%                   default == 0.05
%                   - 'HorizontalAlignment': Horizontal alignment of the text
%                   must be an unambiguous, case-insensitive match to one of:
%                       'left'      - Aligns the left edge of the text with the position
%                       'center'    - Centers the text at the position
%                       'right'     - Aligns the right edge of the text with the position
%                   default == set based on 'TextLocation'
%                   - 'VerticalAlignment': Vertical alignment of the text
%                   must be an unambiguous, case-insensitive match to one of:
%                       'top'       - Aligns the top edge of the text with the position
%                       'cap'       - Aligns with characters that have capital letters
%                       'middle'    - Centers the text vertically at the position
%                       'baseline'  - Aligns with the text baseline
%                       'bottom'    - Aligns the bottom edge of the text with the position
%                   default == set based on 'TextLocation'
%                   - 'AxesHandle': axes handle to plot on
%                   must be an empty or an axes object handle
%                   default == gca
%                   - Any other parameter-value pair for text()
%
% Requires:
%       cd/set_axes_properties.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_correlation_coefficient.m
%       cd/plot_regression_line.m

% File History:
% 2025-09-18 - Created by extracting code from plot_regression_line.m

%% Hard-coded parameters
validTextLocations = {'topleft', 'topright', 'bottomleft', 'bottomright'};
validUnits = {'normalized', 'data', 'pixels', 'inches', 'centimeters', 'points'};

%% Default values for optional arguments
delimiterDefault = '\n';
unitsDefault = 'normalized';
textLocationDefault = 'topleft';
textMarginDefault = 0.05;
hAlignDefault = '';         % Set later based on textLocation
vAlignDefault = '';         % Set later based on textLocation
axHandleDefault = [];       % gca by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;    % allow extraneous options for text()

% Add required inputs to the Input Parser
addRequired(iP, 'textStr', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['textStr must be a character array, a string array, ', ...
            'or a cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Units', unitsDefault, ...
    @(x) any(validatestring(x, validUnits)));
addParameter(iP, 'TextLocation', textLocationDefault, ...
    @(x) any(validatestring(x, validTextLocations)));
addParameter(iP, 'TextMargin', textMarginDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'HorizontalAlignment', hAlignDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'VerticalAlignment', vAlignDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, textStr, varargin{:});
delimiter = iP.Results.Delimiter;
units = validatestring(iP.Results.Units, validUnits);
textLocation = validatestring(iP.Results.TextLocation, validTextLocations);
textMargin = iP.Results.TextMargin;
hAlign = iP.Results.HorizontalAlignment;
vAlign = iP.Results.VerticalAlignment;
axHandle = iP.Results.AxesHandle;

% Keep unmatched arguments for the text() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the axes to plot on
axHandle = set_axes_properties('AxesHandle', axHandle);

% If textStr is a cell array, join it into a single character vector
if iscell(textStr)
    textStr = strjoin(textStr, delimiter);
end

% Don't plot if the string is empty
if isempty(textStr)
    htext = gobjects;
    return;
end

% Get axes limits based on the selected units
if strcmpi(units, 'normalized')
    xLim = [0, 1];
    yLim = [0, 1];

    % For normalized units, margin is a fraction of the range
    xMargin = textMargin * diff(xLim);
    yMargin = textMargin * diff(yLim);
else
    xLim = get(axHandle, 'XLim');
    yLim = get(axHandle, 'YLim');

    % For data units, margin is an absolute value
    xMargin = textMargin;
    yMargin = textMargin;
end

% Determine text position based on limits and margin
switch textLocation
    case 'topleft'
        xPos = xLim(1) + xMargin;
        yPos = yLim(2) - yMargin;
    case 'topright'
        xPos = xLim(2) - xMargin;
        yPos = yLim(2) - yMargin;
    case 'bottomleft'
        xPos = xLim(1) + xMargin;
        yPos = yLim(1) + yMargin;
    case 'bottomright'
        xPos = xLim(2) - xMargin;
        yPos = yLim(1) + yMargin;
end

% Determine alignment based on location (if not provided by user)
switch textLocation
    case 'topleft'
        if isempty(hAlign); hAlign = 'left'; end
        if isempty(vAlign); vAlign = 'top'; end
    case 'topright'
        if isempty(hAlign); hAlign = 'right'; end
        if isempty(vAlign); vAlign = 'top'; end
    case 'bottomleft'
        if isempty(hAlign); hAlign = 'left'; end
        if isempty(vAlign); vAlign = 'bottom'; end
    case 'bottomright'
        if isempty(hAlign); hAlign = 'right'; end
        if isempty(vAlign); vAlign = 'bottom'; end
end

%% Do the job
% Add the text object to the plot and store its handle
htext = text(axHandle, xPos, yPos, textStr, 'Units', units, ...
             'HorizontalAlignment', hAlign, 'VerticalAlignment', vAlign, ...
             otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%