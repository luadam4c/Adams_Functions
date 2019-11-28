function m3ha_compute_and_plot_geometry2(outFolder, suffixes, diamSoma, LDend1, LDend2, linewidth)
%% Plot the geometry of the cell
% TODO: Input parser
% TODO: Make sure that if suffixes is a cell array
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
addRequired(iP, 'suffixes', ...         % corresponding suffixes
	@(x) validateattributes(x, {'char', 'string', 'cell'}, {'nonempty'}));

parse(iP, outFolder, suffixes);
outFolder = iP.Results.outFolder;
suffixes = iP.Results.suffixes;
validfn = @(x) validateattributes(x, {'numeric'}, ...
                                    {'vector', 'numel', numel(suffixes)});
addRequired(iP, 'diamSoma', validfn);
addRequired(iP, 'LDend1', validfn);
addRequired(iP, 'LDend2', validfn);
parse(iP, outFolder, suffixes, diamSoma, LDend1, LDend2);
diamSoma = iP.Results.diamSoma;
Ldend1 = iP.Results.LDend1;
Ldend2 = iP.Results.LDend2;

% Deal with arguments TODO
if iscell(suffixes)         % many conditions are passed in
    % Concatenate all suffixes together into a new suffix
    nSuffixes = numel(suffixes);
    suffix = '';
    for iSuffix = 1:nSuffixes
        suffix = [suffix, suffixes{iSuffix}];
    end
else                        % only one condition is passed in
    nSuffixes = 1;
    suffix = suffixes;
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
for iSuffix = 1:nSuffixes
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
colorTemp = zeros(nSuffixes, 1);
suffixesClean = cell(nSuffixes, 1);
for iSuffix = 1:nSuffixes
	colorTemp(iSuffix) = plot(NaN, NaN, 'LineWidth', linewidth, ...
	                            'Color', color(iSuffix,:));
    suffixesClean{iSuffix} = strrep(suffixes{iSuffix}, '_', '');
end
legend(colorTemp, suffixesClean, 'location', 'eastoutside');
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

    % Make suffixes into a cell array
    suffixes = cell;
    suffixes{1} = suffix;
    @(x) assert((iscell(x) && (min(cellfun(@ischar, x)) ...
                                || min(cellfun(@isstring, x)))) ||
				(isscalar(x) && ischar(x)), ...
        'suffixes must be a cell array of strings/character arrays or scalar text!'));
legend(suffixes, 'location', 'northwest');

hfig.Visible = 'off';

diamDend1 = diamSoma .* diamDendToSoma;
diamDend2 = diamSoma .* diamDendToSoma;
LDend1 = LDend .* (1 - distDendPercent./100);
LDend2 = LDend .* distDendPercent./100;

outFolder = iP.Results.outFolder;
suffixes = iP.Results.suffixes;
validfn = @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', numel(suffixes)});
addRequired(iP, 'diamSoma', validfn);
addRequired(iP, 'LDend', validfn);
addRequired(iP, 'distDendPercent', validfn);
parse(iP, outFolder, suffixes, diamSoma, LDend, diamDendToSoma, distDendPercent);
diamSoma = iP.Results.diamSoma;
Ldend = iP.Results.LDend;
diamDendToSoma = iP.Results.diamDendToSoma;
distDendPercent = iP.Results.distDendPercent;

%}
