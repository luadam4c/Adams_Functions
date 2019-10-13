function m3ha_compute_and_plot_geometry2(outFolder, suffices, diamSoma, LDend1, LDend2, linewidth)
%% Plot the geometry of the cell
% TODO: Input parser
% TODO: Make sure that if suffices is a cell array
%       then diamSoma, LDend, diamDendToSoma, distDendPercent are also arrays of the same length

% File History:
% 2017-09-01 BT - Implemented Input Parser
% 2018-03-08 AL - Made figures invisible upon creation
% 2018-08-15 AL - Updated geometry parameters

color = [0 0 1	% blue
		 1 0 0	% red
		 0 1 0	% green
		 1 0 1];	% purple

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'outFolder', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(iP, 'suffices', ...         % corresponding suffices
	@(x) validateattributes(x, {'char', 'string', 'cell'}, {'nonempty'}));

parse(iP, outFolder, suffices);
outFolder = iP.Results.outFolder;
suffices = iP.Results.suffices;
validfn = @(x) validateattributes(x, {'numeric'}, ...
                                    {'vector', 'numel', numel(suffices)});
addRequired(iP, 'diamSoma', validfn);
addRequired(iP, 'LDend1', validfn);
addRequired(iP, 'LDend2', validfn);
parse(iP, outFolder, suffices, diamSoma, LDend1, LDend2);
diamSoma = iP.Results.diamSoma;
Ldend1 = iP.Results.LDend1;
Ldend2 = iP.Results.LDend2;

% Deal with arguments TODO
if iscell(suffices)         % many conditions are passed in
    % Concatenate all suffices together into a new suffix
    nSuffices = numel(suffices);
    suffix = '';
    for iSuffix = 1:nSuffices
        suffix = [suffix, suffices{iSuffix}];
    end
else                        % only one condition is passed in
    nSuffices = 1;
    suffix = suffices;
end

% Create figure
hfig = figure('Visible', 'off');
clf(hfig);
hold on;
fprintf('Plotting geometry for %s ... \n\n', suffix(2:end));

% Use equal data unit lengths along each axis
daspect([1 1 1]);

% Compute geometry parameters
diamDend1 = diamSoma .* 0.5;
diamDend2 = diamSoma .* 0.3;

% Plot geometry
for iSuffix = 1:nSuffices
    % Plot soma
    rectangle('Position', [-diamSoma(iSuffix)/2, -diamSoma(iSuffix)/2, ...
                            diamSoma(iSuffix), diamSoma(iSuffix)], ...
		     'Curvature', [0, 0], 'EdgeColor', color(iSuffix,:), ...
             'LineWidth', linewidth);

    % Plot dendrite 1
    rectangle('Position', [diamSoma(iSuffix)/2, -diamDend1(iSuffix)/2, ...
                            LDend1(iSuffix), diamDend1(iSuffix)], ...
		     'Curvature', [0, 0], 'EdgeColor', color(iSuffix,:), ...
             'LineWidth', linewidth);

    % Plot dendrite 2
    rectangle('Position', [diamSoma(iSuffix)/2 + LDend1(iSuffix), ...
                            -diamDend2(iSuffix)/2, ...
                            LDend2(iSuffix), diamDend2(iSuffix)], ...
		     'Curvature', [0, 0], 'EdgeColor', color(iSuffix,:), ...
             'LineWidth', linewidth);
end

% Create title, legend and save figure
title(['Geometry for ', strrep(suffix(2:end), '_', '\_')]);
hold on;
colorTemp = zeros(nSuffices, 1);
sufficesClean = cell(nSuffices, 1);
for iSuffix = 1:nSuffices
	colorTemp(iSuffix) = plot(NaN, NaN, 'LineWidth', linewidth, ...
	                            'Color', color(iSuffix,:));
    sufficesClean{iSuffix} = strrep(suffices{iSuffix}, '_', '');
end
legend(colorTemp, sufficesClean, 'location', 'eastoutside');
% set(gca, 'FontSize', 20);
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

    % Make suffices into a cell array
    suffices = cell;
    suffices{1} = suffix;
    @(x) assert((iscell(x) && (min(cellfun(@ischar, x)) ...
                                || min(cellfun(@isstring, x)))) ||
				(isscalar(x) && ischar(x)), ...
        'suffices must be a cell array of strings/character arrays or scalar text!'));
legend(suffices, 'location', 'northwest');

hfig.Visible = 'off';

diamDend1 = diamSoma .* diamDendToSoma;
diamDend2 = diamSoma .* diamDendToSoma;
LDend1 = LDend .* (1 - distDendPercent./100);
LDend2 = LDend .* distDendPercent./100;

outFolder = iP.Results.outFolder;
suffices = iP.Results.suffices;
validfn = @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', numel(suffices)});
addRequired(iP, 'diamSoma', validfn);
addRequired(iP, 'LDend', validfn);
addRequired(iP, 'distDendPercent', validfn);
parse(iP, outFolder, suffices, diamSoma, LDend, diamDendToSoma, distDendPercent);
diamSoma = iP.Results.diamSoma;
Ldend = iP.Results.LDend;
diamDendToSoma = iP.Results.diamDendToSoma;
distDendPercent = iP.Results.distDendPercent;

%}
