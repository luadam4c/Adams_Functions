function [h] = plot_traces(tVec, data, xlimits, xLabel, yLabel, ...
                            traceLabels, figTitle, figName, figNum)
%% Plots traces all in one place
% Usage: [h] = plot_traces(tVec, data, xlimits, xLabel, yLabel, ...
%                           traceLabels, figTitle, figName, figNum)
% Outputs:
%       h           - figure handle for the created figure
%                   must be a figure handle
%
% Arguments: TODO
%       tVec        - time vector for plotting
%       data        - data array (each column is a data vector)
%       xlimits     - x-axis limits
%       xLabel      - x-axis label
%       yLabel      - y-axis label
%       traceLabels - legend labels for each trace
%       figTitle    - figure title
%       figName     - figure name
%       figNum      - figure number
%
% Used by:
%       cd/plot_traces_abf.m

% File History:
% 2018-09-18 Moved from plot_traces_abf.m
% TODO: Combine with plot_signals.m by defining a parameter PlotMode


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of traces to plot
nTraces = size(data, 2);

% Check if traceLabels has the correct length
if numel(traceLabels) ~= nTraces
    error('traceLabels has %d elements instead of %d!!', ...
            numel(traceLabels), nTraces);
end

% Compute minimim and maximum Y values
minY = min(min(data));
maxY = max(max(data));
rangeY = maxY - minY;

% Create a figure and clear it
% h = figure(figNum);
% set(h, 'Visible', 'Off');
h = figure('Visible', 'Off');
clf(h);

% Plot all traces
hold on;
for iTrace = 1:nTraces
    plot(tVec, data(:, iTrace), ...
        'DisplayName', traceLabels{iTrace});
end

% Determine the appropriate axis limits
xlim(xlimits);
if rangeY ~= 0
    ylim([minY - 0.2 * rangeY, maxY + 0.2 * rangeY]);
end

% Generate an x-axis label
xlabel(xLabel);

% Generate a y-axis label
ylabel(yLabel);

% Generate a title
title(figTitle, 'Interpreter', 'none');

% Generate a legend if there is more than one trace
if nTraces > 1 && nTraces < 10
    legend('location', 'northeast');
elseif nTraces >= 10
    legend('location', 'eastoutside');
end

% Save figure
saveas(h, figName, 'png');

% Hold off and close figure
hold off;
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
