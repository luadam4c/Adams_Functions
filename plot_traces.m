function [fig, subPlots, plotsData, plotsDataToCompare] = ...
                plot_traces (tVecs, data, varargin)
%% Plots traces all in one place, overlapped or in parallel
% Usage: [fig, subPlots, plotsData, plotsDataToCompare] = ...
%               plot_traces (tVecs, data, varargin)
% Outputs:
%       fig         - figure handle for the created figure
%                   specified as a figure object handle
%       subPlots    - axes handles for the subplots
%                   specified as a vector of axes object handles
%       plotsData   - line handles for the data plots
%                   specified as a vector of chart line object handles
%       plotsDataToCompare  - line handles for the data to compare plots
%                   specified as a vector of chart line object handles
%
% Arguments:
%       tVecs       - time vector(s) for plotting
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       data        - data vectors(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'OverWrite': whether to overwrite existing output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotMode': plotting mode for multiple traces
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'overlapped'    - overlapped in a single plot
%                       'parallel'      - in parallel in subPlots
%                   must be consistent with plot_traces_abf.m
%                   default == 'overlapped'
%                   - 'DataToCompare': data vector(s) to compare against
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%                   default == []
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_xlimits.m
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_ylimits.m
%                   - 'LinkAxesOption': option for the linkaxes()
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'none' - don't apply the function
%                       'x'    - link x axes only
%                       'y'    - link y axes only
%                       'xy'   - link x and y axes
%                       'off'  - unlink axes
%                   must be consistent with linkaxes()
%                   default == 'none'
%                   - 'XUnits': x-axis units
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'unit'
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == ['Time (', xUnits, ')']
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'Data' if plotMode is 'overlapped'
%                               {'Trace #1', 'Trace #2', ...}
%                                   if plotMode is 'parallel'
%                   - 'TraceLabels': labels for the traces, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Trace #1', 'Trace #2', ...}
%                   - 'ColorMap': a color map that also groups traces
%                                   each set of traces will be on the same row
%                                   if plot mode is 'parallel'
%                   must be a numeric array with 3 columns
%                   default == colormap(jet(nTraces))
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress' if nTraces == 1 
%                               'northeast' if nTraces is 2~9
%                               'eastoutside' if nTraces is 10+
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Traces for ', figName]
%                               or [yLabel, ' over time']
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       cd/apply_iteratively.m
%       cd/argfun.m
%       cd/compute_xlimits.m
%       cd/compute_ylimits.m
%       cd/count_vectors.m
%       cd/create_colormap.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/isfigtype.m
%       cd/islegendlocation.m
%       cd/ispositiveintegerscalar.m
%       cd/match_format_vector_sets.m
%       cd/save_all_figtypes.m
%       ~/Downloaded_Function/suplabel.m
%       ~/Downloaded_Function/subplotsqueeze.m
%
% Used by:
%       cd/m3ha_plot_individual_traces.m
%       cd/plot_traces_abf.m

% File History:
% 2018-09-18 Moved from plot_traces_abf.m
% 2018-09-25 Implemented the input parser
% 2018-09-25 Added 'PlotMode' and 'LegendLocation' as optional parameters
% 2018-10-29 Added 'ColorMap' as an optional parameter
% 2018-10-29 Number of rows in parallel mode is now dependent on the 
%               number of rows in the colorMap provided
% 2018-10-29 Added 'DataToCompare' as an optional parameter
% 2018-10-31 Now uses match_format_vector_sets.m
% 2018-11-01 Now returns axes handles for subplots
% 2018-11-22 Now accepts xLimits as a cell array
% 2018-11-22 Added 'XUnits' as an optional parameter
% 2018-12-15 Fixed the passing of parameters to the helper function
% 2018-12-15 Now returns the axes handle as the second output
%               for overlapped plots
% 2018-12-17 Now uses create_labels_from_numbers.m
% 2018-12-17 Now uses iP.Unmatched
% 2018-12-17 Now uses compute_xlimits.m and compute_ylimits.m
% 2018-12-19 Now returns line object handles for the plots
% 2018-12-19 Added 'FigHandle' as an optional argument
% 2018-12-19 Now restricts vectors to x limits first

%% Hard-coded parameters
validPlotModes = {'overlapped', 'parallel'};
validLinkAxesOptions = {'none', 'x', 'y', 'xy', 'off'};
maxRowsWithOneOnly = 8;
maxNTracesForAnnotations = 8;
maxNTracesForLegends = 12;
subPlotSqeezeFactor = 1.2;

%% Default values for optional arguments
verboseDefault = true;
overWriteDefault = true;        % overwrite previous plots by default
plotModeDefault = 'overlapped'; % plot traces overlapped by default
dataToCompareDefault = [];      % no data to compare against by default
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later
linkAxesOptionDefault = 'none'; % don't force link axes by default
xUnitsDefault = 'unit';         % the default x-axis units
xLabelDefault = '';             % set later
yLabelDefault = '';             % set later
colorMapDefault = [];           % set later
traceLabelsDefault = '';        % set later
legendLocationDefault = 'auto'; % set later
figTitleDefault = '';           % set later
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figNameDefault = '';            % don't save figure by default
figTypesDefault = 'png';        % save as png file by default

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

% Add required inputs to an Input Parser
addRequired(iP, 'tVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'data', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OverWrite', overWriteDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMode', plotModeDefault, ...
    @(x) any(validatestring(x, validPlotModes)));
addParameter(iP, 'DataToCompare', dataToCompareDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || iscell(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'LinkAxesOption', linkAxesOptionDefault, ...
    @(x) any(validatestring(x, validLinkAxesOptions)));
addParameter(iP, 'XUnits', xUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'TraceLabels', traceLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault, ...
    @(x) isempty(x) || isnumeric(x) && size(x, 2) == 3);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, tVecs, data, varargin{:});
verbose = iP.Results.Verbose;
overWrite = iP.Results.OverWrite;
plotMode = validatestring(iP.Results.PlotMode, validPlotModes);
dataToCompare = iP.Results.DataToCompare;
xUnits = iP.Results.XUnits;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
linkAxesOption = validatestring(iP.Results.LinkAxesOption, ...
                                validLinkAxesOptions);
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
traceLabels = iP.Results.TraceLabels;
colorMap = iP.Results.ColorMap;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

%% Preparation
% If data is empty, return
if isempty(data) || iscell(data) && all(cellfun(@isempty, data))
    fprintf('Nothing to plot!\n');
    return
end

% If not to overwrite, check if the figure already exists
if ~overWrite && check_fullpath(figName, 'Verbose', verbose)
    % Skip this figure
    fprintf('%s skipped!\n', figName);
    return;
end

% Restrict to x limits for faster processing
if ~isempty(xLimits) && isnumeric(xLimits)
    % Find the end points
    endPoints = find_window_endpoints(xLimits, tVecs);

    % Restrict to these end points
    [tVecs, data, dataToCompare] = ...
        argfun(@(x) extract_subvectors(x, 'EndPoints', endPoints), ...
                tVecs, data, dataToCompare);
end

% Match the number of vectors between data and dataToCompare
%   and make sure boths are column cell arrays of column vectors
[data, dataToCompare] = ...
    match_format_vector_sets(data, dataToCompare, 'ForceCellOutputs', true);

% Extract number of traces
nTraces = count_vectors(data);

% Decide on the colormap
if isempty(colorMap)
    if nTraces <= maxRowsWithOneOnly
        colorMap = create_colormap(nTraces);
    else
        colorMap = create_colormap(floor(sqrt(nTraces)));
    end
end

% Determine the number of rows and the number of traces per row
nRows = size(colorMap, 1);
nTracesPerRow = ceil(nTraces / nRows);

% Compute minimum and maximum Y values
% TODO: Consider dataToCompare range too
minY = apply_iteratively(@min, [data; dataToCompare]);
maxY = apply_iteratively(@max, [data; dataToCompare]);
rangeY = maxY - minY;

% Force as column cell array and match up to nTraces elements 
tVecs = match_format_vector_sets(tVecs, data);

% Set the default time axis limits
if isempty(xLimits)
    % TODO: Use extract_element.m 'first' for cell arrays?
    if iscell(tVecs)
        tVec = tVecs{1};
    else
        tVec = tVecs;
    end

    xLimits = compute_xlimits(tVec, 'Coverage', 100);
end

% Set the default y-axis limits
if isempty(yLimits) && ~strcmpi(plotMode, 'parallel') && rangeY ~= 0
    % TODO: Deal with yLimits if it is a cell array
    yLimits = compute_ylimits(minY, maxY, 'Coverage', 80);
end

% Set the default x-axis labels
if isempty(xLabel)
    xLabel = ['Time (', xUnits, ')'];
end

% Set the default y-axis labels
if isempty(yLabel)
    switch plotMode
    case 'overlapped'
        yLabel = 'Data';
    case 'parallel'
        if nTraces > 1
            yLabel = create_labels_from_numbers(1:nTraces, 'Prefix', 'Trace #');
        else
            yLabel = {'Data'};
        end
    otherwise
        error(['The plot mode ', plotMode, ' has not been implemented yet!']);
    end
end

% Make sure y-axis labels are consistent
switch plotMode
case 'overlapped'
    if iscell(yLabel)
        fprintf('Only the first yLabel will be used!\n');
        yLabel = yLabel{1};
    end
case 'parallel'
    % Force as column cell array and match up to nTraces elements
    yLabel = match_format_vector_sets(yLabel, data);
otherwise
    error(['The plot mode ', plotMode, ' has not been implemented yet!']);
end

% Set the default trace labels
if isempty(traceLabels)
    traceLabels = create_labels_from_numbers(1:nTraces, 'Prefix', 'Trace #');
end

% Make sure trace labels are cell arrays
if ~isempty(traceLabels) && ...
    (ischar(traceLabels) || isstring(traceLabels)) && ...
    ~strcmpi(traceLabels, 'suppress')
    traceLabels = {traceLabels};
end

% Check if traceLabels has the correct length
if iscell(traceLabels) && numel(traceLabels) ~= nTraces
    error('traceLabels has %d elements instead of %d!!', ...
            numel(traceLabels), nTraces);
end

% Set the default figure title
if isempty(figTitle)
    if ~isempty(figName) && nTraces == 1
        figTitle = ['Traces for ', traceLabels{1}];
    elseif ~isempty(figName)
        figTitle = ['Traces for ', figName];
    elseif ischar(yLabel)
        figTitle = [yLabel, ' over ', xLabel];
    else
        figTitle = ['Data over ', xLabel];        
    end
end

% Set legend location based on number of traces
if strcmpi(legendLocation, 'auto')
    if nTraces > 1 && nTraces <= maxNTracesForAnnotations
        legendLocation = 'northeast';
    elseif nTraces > maxNTracesForAnnotations && nTraces <= maxNTracesForLegends
        legendLocation = 'eastoutside';
    else
        legendLocation = 'suppress';
    end
end

%% Plot data over all possible intervals
if iscell(xLimits)
    % Count the number of intervals
    nIntervals = numel(xLimits);

    % Run through all intervals
    parfor iInterval = 1:nIntervals
    % Get the current x-axis limits
        xLimitsThis = xLimits{iInterval};

        % Create a string for the interval
        intervalStrThis = sprintf('%.0f-%.0f%s', ...
                            xLimitsThis(1), xLimitsThis(2), xUnits);

        % Construct a file suffix
        suffixThis = sprintf('_%s.png', intervalStrThis);

        % Create a new figure name
        figNameThis = replace(figName, '.png', suffixThis);

        % If not to overwrite, check if the figure already exists
        if ~overWrite && check_fullpath(figNameThis, 'Verbose', verbose)
            % Skip this figure
            fprintf('%s skipped!\n', figNameThis);
        else
            % Print to standard output
            if verbose
                fprintf('Interval to show = %s\n', intervalStrThis);
            end
            
            % Create a new figure title
            if ~strcmpi(figTitle, 'suppress')
                figTitleThis = [figTitle, ' (', intervalStrThis, ')'];
            else
                figTitleThis = 'suppress';
            end

            % Find the corresponding index endpoints
            endPoints = find_window_endpoints(xLimitsThis, tVecs, ...
                                                'BoundaryMode', 'inclusive');

            % Truncate all traces
            [tVecsThis, dataThis, dataToCompareThis] = ...
                argfun(@(x) extract_subvectors(x, 'EndPoints', endPoints), ...
                        tVecs, data, dataToCompare);

            % Plot all traces
            fig = plot_traces_helper(verbose, plotMode, ...
                            tVecsThis, dataThis, dataToCompareThis, ...
                            xUnits, xLimitsThis, yLimits, linkAxesOption, ...
                            xLabel, yLabel, traceLabels, colorMap, ...
                            legendLocation, figTitleThis, ...
                            figHandle, figNumber, figNameThis, figTypes, ...
                            nTraces, nRows, nTracesPerRow, ...
                            maxNTracesForAnnotations, subPlotSqeezeFactor, ...
                            otherArguments);
            
            % Hold off and close figure
            hold off;
            close(fig)
        end
    end

    % Return nothing
    fig = gobjects(1);
    subPlots = gobjects(1);
    plotsData = gobjects(1);
    plotsDataToCompare = gobjects(1);
else
    % Plot all traces
    [fig, subPlots, plotsData, plotsDataToCompare] = ...
        plot_traces_helper(verbose, plotMode, ...
                        tVecs, data, dataToCompare, ...
                        xUnits, xLimits, yLimits, linkAxesOption, ...
                        xLabel, yLabel, traceLabels, colorMap, ...
                        legendLocation, figTitle, ...
                        figHandle, figNumber, figName, figTypes, ...
                        nTraces, nRows, nTracesPerRow, ...
                        maxNTracesForAnnotations, subPlotSqeezeFactor, ...
                        otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig, subPlots, plotsData, plotsDataToCompare] = ...
                plot_traces_helper (verbose, plotMode, ...
                        tVecs, data, dataToCompare, ...
                        xUnits, xLimits, yLimits, linkAxesOption, ...
                        xLabel, yLabel, traceLabels, colorMap, ...
                        legendLocation, figTitle, ...
                        figHandle, figNumber, figName, figTypes, ...
                        nTraces, nRows, nTracesPerRow, ...
                        maxNTracesForAnnotations, subPlotSqeezeFactor, ...
                        otherArguments)

% Decide on the figure to plot on
if ~isempty(figHandle)
    fig = figure(figHandle);
elseif ~isempty(figNumber)
    fig = figure(figNumber);
else
    fig = gcf;
end

% Clear the figure
clf(fig);

% Initialize graphics object arrays for plots
plotsData = gobjects(nTraces, 1);
plotsDataToCompare = gobjects(nTraces, 1);

switch plotMode
case 'overlapped'
    % Hold on
    hold on

    % Plot all traces together
    for iTrace = 1:nTraces
        % Get the current row (color) number
        thisRowNumber = ceil(iTrace/nTracesPerRow);

        % Plot data to compare against as a black trace
        if ~isempty(dataToCompare{iTrace})
            plotsDataToCompare(iTrace) = ...
                plot(tVecs{iTrace}, dataToCompare{iTrace}, ...
                                'Color', 'k', otherArguments);
        end
        
        % Plot the data using the color map
        p = plot(tVecs{iTrace}, data{iTrace}, ...
                'Color', colorMap(thisRowNumber, :), otherArguments);

        % Set the legend label as the trace label if provided
        if ~strcmpi(traceLabels, 'suppress')
            set(p, 'DisplayName', traceLabels{iTrace});
        end

        % Store handles in array
        plotsData(iTrace) = p;
    end
    
    % Set time axis limits
    if ~iscell(xLimits) && ~strcmpi(xLimits, 'suppress')
        xlim(xLimits);
    end

    % Set y axis limits
    if ~isempty(yLimits) && ~strcmpi(yLimits, 'suppress')
        ylim(yLimits);
    end

    % Generate an x-axis label
    if ~strcmpi(xLabel, 'suppress')
        xlabel(xLabel);
    end

    % Generate a y-axis label
    if ~strcmpi(yLabel, 'suppress')
        ylabel(yLabel);
    end

    % Generate a title
    if ~strcmpi(figTitle, 'suppress')
        title(figTitle, 'Interpreter', 'none');
    end

    % Generate a legend if there is more than one trace
    if ~strcmpi(legendLocation, 'suppress')
        legend(p1, 'location', legendLocation);
    end

    % Save current axes handle
    subPlots = gca;
case 'parallel'
    if ~strcmpi(legendLocation, 'suppress')
        % Set a legend location differently    
        legendLocation = 'northeast';
    end

    % Initialize graphics object arrays for subplots
    subPlots = gobjects(nTraces, 1);

    % Plot each trace as a different subplot
    %   Note: the number of rows is based on the number of rows in the color map
    for iTrace = 1:nTraces
        % Create a subplot and hold on
        ax = subplot(nRows, nTracesPerRow, iTrace); hold on

        % Get the current row number
        thisRowNumber = ceil(iTrace/nTracesPerRow);

        % Get the current column number
        if nTracesPerRow > 1
            thisColNumber = mod(iTrace, nTracesPerRow);
        else
            thisColNumber = 1;
        end
        
        % Plot data to compare against as a black trace
        if ~isempty(dataToCompare{iTrace})
            plotsDataToCompare(iTrace) = ...
                plot(tVecs{iTrace}, dataToCompare{iTrace}, ...
                        'Color', 'k', otherArguments);
        end

        % Plot the data using the color map
        p = plot(tVecs{iTrace}, data{iTrace}, ...
                    'Color', colorMap(thisRowNumber, :), otherArguments);

        % Set the legend label as the trace label if provided
        if ~strcmpi(traceLabels, 'suppress')
            set(p, 'DisplayName', traceLabels{iTrace});
        end

        % Set time axis limits
        if ~iscell(xLimits) && ~strcmpi(xLimits, 'suppress')
            xlim(xLimits);
        end

        % Set y axis limits
        if ~isempty(yLimits) && ~strcmpi(yLimits, 'suppress')
            ylim(yLimits);
        end

        % Generate a y-axis label
        if ~strcmpi(yLabel{iTrace}, 'suppress')
            ylabel(yLabel{iTrace});
        end

        % Generate a legend
        if ~strcmpi(legendLocation, 'suppress')
            legend(p1, 'location', legendLocation);
        end

        % Remove x tick labels except for the last row
        if thisRowNumber ~= nRows
            set(ax, 'XTickLabel', []);
        end

        % Remove x tick labels except for the first column
        if thisColNumber ~= 1
            set(ax, 'YTickLabel', []);
        end

        % Create a title for the first subplot
        if ~strcmpi(figTitle, 'suppress') && ...
            nTracesPerRow == 1 && iTrace == 1
            title(figTitle, 'Interpreter', 'none');
        end

        % Create a label for the X axis only for the last row
        if ~strcmpi(xLabel, 'suppress') && nTracesPerRow == 1 && ...
            iTrace == nTraces
            xlabel(xLabel);
        end

        % Store handles in array
        subPlots(iTrace) = ax;
        plotsData(iTrace) = p;
    end

    % If requested, link or unlink axes of subPlots
    if ~strcmpi(linkAxesOption, 'none')
        linkaxes(subPlots, linkAxesOption);
    end

    % If nTraces > maxNTracesForAnnotations, expand all subPlots by 1.2
    if nTraces > maxNTracesForAnnotations
        subplotsqueeze(fig, subPlotSqeezeFactor);
    end
    
    % Create an overarching title
    if ~strcmpi(figTitle, 'suppress') && nTracesPerRow > 1
        suptitle(figTitle);
    end

    % Create an overarching x-axis label
    if ~strcmpi(xLabel, 'suppress') && nTracesPerRow > 1
        suplabel(xLabel, 'x');
    end
otherwise
    error(['The plot mode ', plotMode, ' has not been implemented yet!']);
end

%% Save
% Save figure
if ~isempty(figName)
    % TODO: Save figure with other varying attributes
    if iscell(xLimits)
        % TODO: Pull out to function save_all_intervals.m
        %   Note: this part is very slow for large data

        % Count the number of intervals
        nIntervals = numel(xLimits);

        % Run through all intervals
        for iInterval = 1:nIntervals
            % Get the current x-axis limits
            xLimitsThis = xLimits{iInterval};

            % Create a string for the interval
            intervalStrThis = sprintf('%.0f-%.0f%s', ...
                                xLimitsThis(1), xLimitsThis(2), xUnits);

            % Print to standard output
            if verbose
                fprintf('Interval to show = %s\n', intervalStrThis);
            end
            
            % Create a new figure title
            if ~strcmpi(figTitle, 'suppress')
                figTitleThis = [figTitle, ' (', intervalStrThis, ')'];
            else
                figTitleThis = 'suppress';
            end

            % Change the x-axis limits
            switch plotMode
            case 'overlapped'
                % Change the figure title
                if ~strcmpi(figTitleThis, 'suppress')
                    title(figTitleThis, 'Interpreter', 'none');
                end

                % Change the x-axis limits
                xlim(xLimitsThis);
            case 'parallel'
                for iTrace = 1:nTraces
                    % Go to the subplot
                    subplot(subPlots(iTrace));

                    % Create a title for the first subplot
                    if ~strcmpi(figTitleThis, 'suppress') && ...
                        nTracesPerRow == 1 && iTrace == 1
                        title(figTitleThis, 'Interpreter', 'none');
                    end

                    % Change x-axis limits
                    xlim(xLimitsThis);
                end

                % Create an overarching title
                if ~strcmpi(figTitleThis, 'suppress') && nTracesPerRow > 1
                    suptitle(figTitleThis);
                end
            end

            % Construct a file suffix
            suffixThis = sprintf('_%s.png', intervalStrThis);

            % Create a new figure name
            figNameThis = replace(figName, '.png', suffixThis);

            % Save the new figure
            save_all_figtypes(fig, figNameThis, figTypes);
        end
    else
        % Save the new figure
        save_all_figtypes(fig, figName, figTypes);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
OLD CODE:

function fig = plot_traces(tVec, data, xLimits, xLabel, yLabel, ...
                            traceLabels, figTitle, figName, figNum)
%       xLimits     - x-axis limits
%       xLabel      - x-axis label
%       yLabel      - y-axis label
%       traceLabels - legend labels for each trace
%       figTitle    - figure title
%       figName     - figure name
%       figNum      - figure number

% Hold off and close figure
hold off;
close(fig);

saveas(fig, figName, 'png');

fig = figure(figNum);
set(fig, 'Visible', 'Off');

% Determine the appropriate time axis limits
if ~isempty(xLimits)
    if ~strcmpi(xLimits, 'suppress')
        xlim(xLimits);
    end
else
    xlim([min(tVec), max(tVec)]);
end

% Determine the appropriate y axis limits
if ~isempty(yLimits)
    if ~strcmpi(yLimits, 'suppress')
        ylim(yLimits);
    end
else
    if rangeY ~= 0
        ylim([minY - 0.2 * rangeY, maxY + 0.2 * rangeY]);
    end
end

    @(x) any(validatestring(x, validLegendLocations)));
legendLocation = validatestring(iP.Results.LegendLocation, ...
                                validLegendLocations);

subplot(nTraces, 1, iTrace);

if ~iscell(yLabel)
    yLabel = {yLabel};
end
if iscell(yLabel)
    if numel(yLabel) > nTraces
        fprintf('Too many y labels! Only some will be used!\n');
    elseif numel(yLabel) < nTraces
        fprintf('Not enough y labels!!\n');
        return;
    end
end

nTraces = size(data, 2);

xLimits = cellfun(@(x), [min(x), max(x)], tVecs, 'UniformOutput', false);

minY = min(min(data));
maxY = max(max(data));

yLabel = cell(1, nTraces);
parfor iTrace = 1:nTraces
    yLabel{iTrace} = ['Trace #', num2str(iTrace)];
end

traceLabels = cell(1, nTraces);
parfor iTrace = 1:nTraces
    traceLabels{iTrace} = ['Trace #', num2str(iTrace)];
end

% Hold on if more than one trace
if nTraces > 1
    hold on
end
% Hold off if more than one trace
if nTraces > 1
    hold off
end

% Force as column cell arrays
yLabel = force_column_cell(yLabel);

% Match up to nTraces elements
yLabel = match_dimensions(yLabel, [nTraces, 1]);

% Force data vectors as column cell arrays of column vectors
[tVecs, data, dataToCompare] = ...
    argfun(@force_column_cell, tVecs, data, dataToCompare);

% Match the number of vectors between data and dataToCompare
[data, dataToCompare] = match_array_counts(data, dataToCompare);

% Match the dimensions of tVecs to data
tVecs = match_dimensions(tVecs, size(data));

%       cd/argfun.m
%       cd/force_column_cell.m
%       cd/match_dimensions.m
%       cd/match_array_counts.m

% Hold off
hold off

% Hold off
hold off

if ~strcmpi(xLimitsThis, 'suppress')
end

% Create a new figure number
if ~isempty(figNumber)
    figNumberThis = figNumber + rand(1) * 1000;
else
    figNumberThis = [];
end

yLabel = arrayfun(@(x) ['Trace #', num2str(x)], ...
                    transpose(1:nTraces), 'UniformOutput', false);
traceLabels = arrayfun(@(x) ['Trace #', num2str(x)], ...
                        transpose(1:nTraces), 'UniformOutput', false);

% Compute minimum and maximum time values
minT = min(cellfun(@min, tVecs));
maxT = max(cellfun(@max, tVecs));

% Compute x limits
xLimits = [minT, maxT];

yLimits = [minY - 0.2 * rangeY, maxY + 0.2 * rangeY];

if ~isempty(figName)
    % Create an invisible figure and clear it
    if ~isempty(figNumber)
        fig = figure(figNumber);
        set(fig, 'Visible', 'off');
    else
        fig = figure('Visible', 'off');
    end
    clf(fig);
else
    % Get the current figure
    fig = gcf;
end

set(fig, 'Visible', 'off');

axes(subPlots(iTrace));

    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

minY = min([cellfun(@min, data), cellfun(@min, dataToCompare)]);
maxY = max([cellfun(@max, data), cellfun(@max, dataToCompare)]);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
