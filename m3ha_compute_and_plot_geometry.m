function m3ha_compute_and_plot_geometry(outFolder, suffix, diamSoma, LDend1, LDend2, color, linewidth)
%% Plot the geometry of the cell
%
% 2018-03-08 Made figures invisible
% 2018-08-15 Updated geometry parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct direct parameters
diamDend1 = diamSoma * 0.5;
diamDend2 = diamSoma * 0.3;

% Create figure
hfig = figure('Visible', 'off');
clf(hfig);
hold on;
fprintf('Plotting geometry for %s ... \n\n', suffix(2:end));

% Plot soma
rectangle('Position', [-diamSoma/2, -diamSoma/2, diamSoma, diamSoma], ...
		'Curvature', [0, 0], 'EdgeColor', color, 'LineWidth', linewidth);

% Plot dendrite 1
rectangle('Position', [diamSoma/2, -diamDend1/2, LDend1, diamDend1], ...
		'Curvature', [0, 0], 'EdgeColor', color, 'LineWidth', linewidth);

% Plot dendrite 2
rectangle('Position', [diamSoma/2 + LDend1, -diamDend2/2, LDend2, diamDend2], ...
		'Curvature', [0, 0], 'EdgeColor', color, 'LineWidth', linewidth);

% Use equal data unit lengths along each axis
daspect([1 1 1]);

% Create title, legend and save figure
title(['Geometry for ', strrep(suffix(2:end), '_', '\_')]);
xlim([-30, 260]);
ylim([-30, 30]);
xlabel('Length (um)');
ylabel('Width (um)');
set(gca, 'FontSize', 20);
saveas(hfig, fullfile(outFolder, ['geometry', suffix]), 'png');
close(hfig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

diamDend1 = diamSoma * diamDendToSoma;
diamDend2 = diamSoma * diamDendToSoma;
LDend1 = LDend * (1 - distDendPercent/100);
LDend2 = LDend * distDendPercent/100;


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%