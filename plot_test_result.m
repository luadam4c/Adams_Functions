function handles = plot_test_result (testPValues, varargin)
%% Plots statistical test results (p-values and significance markers) on a figure
% Usage: handles = plot_test_result (testPValues, varargin)
%
% Explanation:
%       This function overlays statistical test results onto an existing plot.
%       It takes a vector of p-values and parameters for positioning.
%       For each p-value, it displays a text string (e.g., 'p_t = 0.04').
%       The text is colored black for non-significant results, red for
%       significant results, and gray if the test was deemed inappropriate.
%       It also plots a symbol ('*' by default) for significance or 'NS'
%       for non-significance.
%
% Example(s):
%       % Create a base plot
%       figure;
%       plot(1:10, rand(1, 10));
%       xlim([0, 11]); ylim([0, 1.5]);
%
%       % Define some sample test results
%       pvals = [0.04, 0.6, 0.01];
%       x_axis_values = 1:4;
%
%       % Provide a cell array of test functions to generate labels automatically
%       plot_test_result(pvals, 'TestFunction', {'ttest', 'ranksum', 'anova'}, ...
%                        'XValuesAll', x_axis_values);
%
%       % Provide a custom symbol for significance
%       plot_test_result([0.02], 'Symbol', 'd', 'XLocText', [9]);
%
%       % Provide a mixed cell array of p-strings and symbols
%       plot_test_result(pvals, 'PString', {'p_t', 'p-value', 'p_ANOVA'}, ...
%                        'Symbol', {'s', 's', 'p'}, ...
%                        'XLocText', [6, 7, 8], 'XLocStar', [6, 7, 8]);
%
% Outputs:
%       handles     - A structure containing handles to the plotted objects:
%                       .pText - handles to the p-value text objects
%                       .sigMarker - handles to the significance markers
%                   specified as a structure
%
% Arguments:
%       testPValues - p-values to plot
%                   must be a numeric vector
%       varargin    - 'PString': Prefix for the p-value text.
%                   must be a char array, a string, or a cell array of strings
%                   of the same length as testPValues.
%                   default == 'p'
%                   - 'TestFunction': Name(s) of the statistical test (e.g., 'ttest')
%                   must be a char array, a string, or a cell array of strings
%                   of the same length as testPValues.
%                   default == ''
%                   - 'Symbol': Marker to use for significant results.
%                   must be a char array, a string, or a cell array of strings
%                   of the same length as testPValues.
%                   default == '*'
%                   - 'XLocText': Absolute x-coordinates for p-value text
%                   must be a numeric vector
%                   default == [] (calculated from XValuesAll, XLocTextRel)
%                   - 'XLocStar': Absolute x-coordinates for markers
%                   must be a numeric vector
%                   default == [] (calculated from XValuesAll, XLocStarRel)
%                   - 'XValuesAll': The unique x-axis values of the original plot
%                   must be a numeric vector
%                   default == 1:numel(testPValues)
%                   - 'YLocTextRel': Relative vertical position for p-value text (0-1)
%                   must be a numeric scalar
%                   default == 0.2
%                   - 'YLocStarRel': Relative vertical position for marker (0-1)
%                   must be a numeric scalar
%                   default == 0.9
%                   - 'XLocTextRel': Relative horizontal position for p-value text (0-1)
%                   must be a numeric scalar
%                   default == 0.5
%                   - 'XLocStarRel': Relative horizontal position for marker (0-1)
%                   must be a numeric scalar
%                   default == 0.5
%                   - 'SigLevel': The significance level (alpha)
%                   must be a numeric scalar
%                   default == 0.05
%                   - 'IsAppropriate': Whether the test was appropriate for the data
%                   must be a logical vector of the same size as testPValues
%                   default == true for all values
%                   - 'AxesHandle': axes handle to plot on
%                   must be a empty or an axes object handle
%                   default == set in set_axes_properties.m
%                   
% Requires:
%       cd/create_error_for_nargin.m
%       cd/hold_off.m
%       cd/hold_on.m
%
% Used by:
%       cd/plot_grouped_jitter.m
%       cd/plot_tuning_curve.m
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt\virt_moore.m

% File History:
% 2025-08-29 Created by Gemini by extracting from plot_tuning_curve.m
% 2025-08-29 Refactored by Gemini to use inputParser and optional x-locations
% 2025-09-11 Added 'TestFunction' optional argument by Gemini
% 2025-09-11 Allowed 'PString' and 'TestFunction' to be vectors by Gemini
% 2025-09-11 Added 'Symbol' optional argument by Gemini
% 2025-09-17 Added 'AxesHandle' as an optional argument
%

%% Default values for optional arguments
pStringDefault = 'p';
testFunctionDefault = '';
symbolDefault = '';
xLocTextDefault = [];
xLocStarDefault = [];
yLocTextRelDefault = 0.2;
yLocStarRelDefault = 0.9;
xLocTextRelDefault = 0.5;
xLocStarRelDefault = 0.5;
xValuesAllDefault = [];         % set later
sigLevelDefault = 0.05;
isAppropriateDefault = [];      % set later
axHandleDefault = [];           % gca by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true; % Allow extraneous options (though none are used)

% Add required inputs to the Input Parser
addRequired(iP, 'testPValues', @isnumeric);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PString', pStringDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'TestFunction', testFunctionDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Symbol', symbolDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'XLocText', xLocTextDefault, @isnumeric);
addParameter(iP, 'XLocStar', xLocStarDefault, @isnumeric);
addParameter(iP, 'YLocTextRel', yLocTextRelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'YLocStarRel', yLocStarRelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'XLocTextRel', xLocTextRelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'XLocStarRel', xLocStarRelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'XValuesAll', xValuesAllDefault, @isnumeric);
addParameter(iP, 'SigLevel', sigLevelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'IsAppropriate', isAppropriateDefault, @islogical);
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, testPValues, varargin{:});
pString = iP.Results.PString;
testFunction = iP.Results.TestFunction;
symbol = iP.Results.Symbol;
xLocText = iP.Results.XLocText;
xLocStar = iP.Results.XLocStar;
yLocTextRel = iP.Results.YLocTextRel;
yLocStarRel = iP.Results.YLocStarRel;
xLocTextRel = iP.Results.XLocTextRel;
xLocStarRel = iP.Results.XLocStarRel;
xValuesAll = iP.Results.XValuesAll;
sigLevel = iP.Results.SigLevel;
isAppropriate = iP.Results.IsAppropriate;
axHandle = iP.Results.AxesHandle;

%% Preparation
% Decide on the axes to plot on
axHandle = set_axes_properties('AxesHandle', axHandle);

% Get the number of p-values to plot
nValues = numel(testPValues);

% Validate vector argument lengths
if iscell(pString) && numel(pString) ~= nValues
    error('PString must be a scalar or have the same number of elements as testPValues');
end
if iscell(testFunction) && numel(testFunction) ~= nValues
    error('TestFunction must be a scalar or have the same number of elements as testPValues');
end
if iscell(symbol) && numel(symbol) ~= nValues
    error('Symbol must be a scalar or have the same number of elements as testPValues');
end

% Set default for isAppropriate if not provided
if isempty(isAppropriate)
    isAppropriate = true(1, nValues);
end

% Decide on the x locations for text and markers
if isempty(xLocText) && isempty(xLocStar)
    % If neither is provided, calculate both from relative positions
    if isempty(xValuesAll)
        xValuesAll = 1:(nValues + 1);
    end
    xInterval = diff(xValuesAll);
    if isempty(xInterval)
        xInterval = 1; % Default interval if only one p-value
    else
        xInterval = xInterval(1); % Assume uniform spacing
    end
    xLocText = xValuesAll(1:end-1) + xInterval * xLocTextRel;
    xLocStar = xValuesAll(1:end-1) + xInterval * xLocStarRel;
elseif isempty(xLocText)
    % If only star locations are provided, use them for text too
    xLocText = xLocStar;
elseif isempty(xLocStar)
    % If only text locations are provided, use them for stars too
    xLocStar = xLocText;
end
% If both xLocText and xLocStar are provided, they are used as is.

% Get current y axis limits to determine absolute y positions
yLimitsNow = get(axHandle, 'YLim');
yRange = yLimitsNow(2) - yLimitsNow(1);

% Decide on the absolute y location for texts and stars
yLocText = yLimitsNow(1) + yRange * yLocTextRel;
yLocStar = yLimitsNow(1) + yRange * yLocStarRel;

% Pre-allocate handle arrays
pTextHandles = gobjects(nValues, 1);
sigMarkerHandles = gobjects(nValues, 1);

%% Do the job
% Hold on
wasHold = hold_on(axHandle);

% Plot texts and markers for each p-value
for iValue = 1:nValues
    % Get the current values for this iteration
    testPValueThis = testPValues(iValue);
    xLocTextThis = xLocText(iValue);
    xLocStarThis = xLocStar(iValue);
    isAppropriateThis = isAppropriate(iValue);

    % Decide on the symbol for this iteration
    if iscell(symbol)
        symbolToUse = symbol{iValue};
    else
        symbolToUse = symbol;
    end

    % If empty symbol passed in, decide based on significance
    if isempty(symbolToUse)
        if testPValueThis < sigLevel
            symbolToUse = '*';
        else
            symbolToUse = 'NS';
        end
    end

    % --- Logic for determining the p-value string for this iteration ---
    pStringToUse = '';
    pStringIsDefault = true;

    % Check if a specific PString was provided for this value
    if iscell(pString)
        pStringToUse = pString{iValue};
        pStringIsDefault = false;
    elseif ~ismember('PString', iP.UsingDefaults)
        pStringToUse = pString;
        pStringIsDefault = false;
    end

    % If PString was default, check for TestFunction
    if pStringIsDefault
        testFunctionToUse = '';
        if iscell(testFunction)
            testFunctionToUse = testFunction{iValue};
        elseif ~isempty(testFunction)
            testFunctionToUse = testFunction;
        end
        
        if ~isempty(testFunctionToUse)
            pStringToUse = ['p_{', testFunctionToUse, '}'];
        else
            pStringToUse = pStringDefault; % Fallback to the absolute default
        end
    end

    % Create the final p-value string with 2 significant digits
    pValueString = [pStringToUse, ' = ', num2str(testPValueThis, 2)];

    % Determine the color for the text based on appropriateness and significance
    if ~isAppropriateThis
        pColor = [0.5, 0.5, 0.5];       % Gray for inappropriate tests
    elseif testPValueThis < sigLevel
        pColor = 'r';                   % Red for significant results
    else
        pColor = 'k';                   % Black for non-significant results
    end

    % Plot the p-value text
    pTextHandles(iValue) = text(axHandle, xLocTextThis, yLocText, ...
                                pValueString, 'Color', pColor, ...
                                'HorizontalAlignment', 'center');

    % Plot the specified symbol or text
    switch symbolToUse
    case {'*', '**', '***'}
        % Use text() with larger, bold font for star symbols
        sigMarkerHandles(iValue) = text(axHandle, xLocStarThis, yLocStar, ...
                                        symbolToUse, 'Color', pColor, ...
                                        'HorizontalAlignment', 'center', ...
                                        'FontSize', 14, 'FontWeight', 'bold');
    otherwise
        % Use text() for 'NS' or other custom symbols
        sigMarkerHandles(iValue) = text(axHandle, xLocStarThis, yLocStar, ...
                                        symbolToUse, 'Color', pColor, ...
                                        'HorizontalAlignment', 'center');
    end
end

% Hold off
hold_off(wasHold, axHandle);

%% Output results
handles.pText = pTextHandles;
handles.sigMarker = sigMarkerHandles;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

