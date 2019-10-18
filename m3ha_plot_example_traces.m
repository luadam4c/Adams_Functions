pharmLabels = {'Control', 'GAT1', 'GAT3', 'Dual'};
fileNames = {'C101210_0004_14', 'C101210_0005_14', 'C101210_0006_14', 'C101210_0007_14'};
fileNamesToCompare = {'C101210_0004_11', 'C101210_0005_11', 'C101210_0006_11', 'C101210_0007_11'};
stimStartWindow = 1000;
responseWindow = [800, 2000];
figTypes = {'png', 'epsc2'};

[data, sweepInfo] = ...
    m3ha_import_raw_traces(fileNames, 'StimStartWindow', stimStartWindow, ...
                            'ResponseWindow', responseWindow);
[dataToCompare, sweepInfo] = ...
    m3ha_import_raw_traces(fileNamesToCompare, 'StimStartWindow', stimStartWindow, ...
                            'ResponseWindow', responseWindow);

[tVecs, vVecs] = extract_columns(data, 1:2);
[tVecsToCompare, vVecsToCompare] = extract_columns(dataToCompare, 1:2);

[fig, subplots] = ...
    plot_traces(tVecs, vVecs, 'DataToCompare', vVecsToCompare, ...
                'PlotMode', 'parallel', 'LineWidth', 2, ...
                'LinkAxesOption', 'xy', 'FigTitle', 'suppress', ...
                'XLabel', 'suppress', 'YLabel', pharmLabels, ...
                'LegendLocation', 'suppress');

for iSubplot = 1:numel(subplots)
    subplot(subplots(iSubplot));
    plot_vertical_line(1000, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
end

update_figure_for_corel(fig);

save_all_figtypes(fig, 'example_traces', figTypes);
