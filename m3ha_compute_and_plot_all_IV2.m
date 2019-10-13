function [v, IInitIT, IInfIT, IInfIh, IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP] = m3ha_compute_and_plot_all_IV2(varargin)
%% Plot I-V curves of all currents together
% Usage: [v, IInitIT, IInfIT, IInfIh, IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP] = m3ha_compute_and_plot_all_IV2(varargin)
%
% Arguments:    
%       v           (opt) voltage vector (x-values of plots)
%                   must be an increasing numeric vector
%                   default: -120:0.5:30 mV
%       varargin    - 'Celsius': temperature of simulation [degC]
%                   must be a numeric scalar
%                   default == 33 degC
%                   - 'Ek': K+ reversal potential [mV]
%                   must be a numeric scalar
%                   default == -100 mV
%                   - 'Ena': Na+ reversal potential [mV]
%                   must be a numeric scalar
%                   default == 88 mV
%                   - 'Eh': reversal potential of Ih [mV]
%                   must be a numeric scalar
%                   default == -43 mV
%                   - 'CaOut': extracellular [Ca++] [mM]
%                   must be a numeric scalar
%                   default == 2 mM
%                   - 'CaIn': intracellular [Ca++] [mM]
%                   must be a numeric scalar
%                   default == 2.4e-4 mM
%                   - 'PcabarIT': maximum Ca++ permeability [cm/s]
%                   must be a numeric scalar
%                   default == 0.2e-3 cm/s
%                   - 'ShiftmIT': depolarizing shift of activation curves [mV]
%                   must be a numeric scalar
%                   default == 1 mV
%                   - 'ShifthIT': depolarizing shift of inactivation curves [mV]
%                   must be a numeric scalar
%                   default == 1 mV
%                   - 'SlopemIT': scaling factor for slope of activation curves
%                   must be a numeric scalar
%                   default == 1
%                   - 'SlopehIT': scaling factor for slope of inactivation curves
%                   must be a numeric scalar
%                   default == 1
%                   - 'GhbarIh': maximum conductance of Ih [S/cm2]
%                   must be a numeric scalar
%                   default == 2.2e-5 S/cm2
%                   - 'ShiftmIh': depolarizing shift of activation curves [mV]
%                   must be a numeric scalar
%                   default == 0 mV
%                   - 'GkbarIKir': maximum conductance of IKir [S/cm2]
%                   must be a numeric scalar
%                   default == 2.0e-5 S/cm2
%                   - 'GkbarIA': maximum conductance of IA [S/cm2]
%                   must be a numeric scalar
%                   default == 5.5e-3 S/cm2
%                   - 'GnabarINaP': maximum conductance of INaP [S/cm2]
%                   must be a numeric scalar
%                   default == 5.5e-6 S/cm2
%                   - 'OutFolder': output directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Suffices': suffices for figure names (includes underscore)
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'PlotInfFlag': whether to plot steady state curves only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotAllFlag': whether to plot all I-V curves
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Used by:    
%
% Requires:
%       cd/m3ha_compute_and_plot_IT.m
%       cd/m3ha_compute_and_plot_Ih.m
%       cd/m3ha_compute_and_plot_IKir.m
%       cd/m3ha_compute_and_plot_IA.m
%       cd/m3ha_compute_and_plot_INaP.m
%
% File History:
% 2017-09-22 BT - Adapted from m3ha_compute_and_plot_INaP2.m
% 2017-10-04 BT - Inputs must be size == numel(suffices), except celsius and eNa
% 2018-03-08 AL - Made figures invisible upon creation

%% Default voltage vector limits and steps
vMinDefault = -120;         % lowest voltage to plot [mV]
vMaxDefault = 30;           % highest voltage to plot [mV]
vStepDefault = 0.05;        % voltage step for plotting [mV]

%% Default parameters
celsiusDefault = 33;        % default temperature of simulation [degC]
eKDefault = -100;           % default K+ reversal potential [mV]
eNaDefault = 88;            % default Na+ reversal potential [mV]
ehDefault = -43;            % default reversal potential of Ih [mV]
caOutDefault = 2;           % default extracellular [Ca++] [mM]
caInDefault = 2.4e-4;       % default intracellular [Ca++] [mM]
pcabarITDefault = 0.2e-3;   % default maximum Ca++ permeability [cm/s]
shiftmITDefault = 1;  % default depolarizing shift of activation curves [mV]
shifthITDefault = 1;  % default depolarizing shift of inactivation curves [mV]
slopemITDefault = 1;  % default scaling factor for slope of activation curves
slopehITDefault = 1;  % default scaling factor for slope of inactivation curves
ghbarIhDefault = 2.2e-5;    % default maximum conductance of Ih [S/cm2]
shiftmIhDefault = 0;    % default depolarizing shift of activation curves [mV]
gkbarIKirDefault = 2.0e-5;  % default maximum conductance of IKir [S/cm2]
gkbarIADefault = 5.5e-3;    % default maximum conductance of IA [S/cm2]
gnabarINaPDefault = 5.5e-6; % default maximum conductance of INaP [S/cm2]
outFolderDefault = pwd;     % default output directory
sufficesDefault = '';         % default suffices for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves only
plotAllFlagDefault = true;  % whether to plot all I-V curves by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_all_IV2';

% Check if Suffices is scalar or cell array
addParameter(iP, 'Suffices', sufficesDefault, ...
	@(x) validateattributes(x, {'char', 'string', 'cell'}, {'nonempty'}));
potSufficesInd = find(strcmp(varargin, 'Suffices'));	% potential index for Suffices string
if ~isempty(potSufficesInd)
	parse(iP, varargin{potSufficesInd}, varargin{potSufficesInd+1});	% parse suffices if exists
else
	parse(iP);
end
suffices = iP.Results.Suffices;
if ~iscell(suffices) & numel(suffices) == 1	% if scalar text, convert to cell
	suffices = {suffices};
end

% Add optional inputs to the Input Parser
addOptional(iP, 'v', [], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}));

% Add parameter-value pairs to the Input Parser
validfn_suff = @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', numel(suffices)});	% size match to numel(suffices)
validfn_single = @(x) validateattributes(x, {'numeric'}, {'scalar'});	% celsius, eK, eNa, eh, caOut, caIn are single
%celsiusDefault = ones(1,numel(suffices)) * celsiusDefault;
%eKDefault = ones(1,numel(suffices)) * eKDefault;
%eNaDefault = ones(1,numel(suffices)) * eNaDefault;
%ehDefault = ones(1,numel(suffices)) * ehDefault;
%caOutDefault = ones(1,numel(suffices)) * caOutDefault;
%caInDefault = ones(1,numel(suffices)) * caInDefault;
pcabarITDefault = ones(1,numel(suffices)) * pcabarITDefault;
shiftmITDefault = ones(1,numel(suffices)) * shiftmITDefault;
shifthITDefault = ones(1,numel(suffices)) * shifthITDefault;
slopemITDefault = ones(1,numel(suffices)) * slopemITDefault;
slopehITDefault = ones(1,numel(suffices)) * slopehITDefault;
ghbarIhDefault = ones(1,numel(suffices)) * ghbarIhDefault;
shiftmIhDefault = ones(1,numel(suffices)) * shiftmIhDefault;
gkbarIKirDefault = ones(1,numel(suffices)) * gkbarIKirDefault;
gkbarIADefault = ones(1,numel(suffices)) * gkbarIADefault;
gnabarINaPDefault = ones(1,numel(suffices)) * gnabarINaPDefault;
addParameter(iP, 'Celsius', celsiusDefault, validfn_single);
addParameter(iP, 'Ek', eKDefault, validfn_single);
addParameter(iP, 'Ena', eNaDefault, validfn_single);
addParameter(iP, 'Eh', ehDefault, validfn_single);
addParameter(iP, 'CaOut', caOutDefault, validfn_single);
addParameter(iP, 'CaIn', caInDefault, validfn_single);
addParameter(iP, 'PcabarIT', pcabarITDefault, validfn_suff);
addParameter(iP, 'ShiftmIT', shiftmITDefault, validfn_suff);
addParameter(iP, 'ShifthIT', shifthITDefault, validfn_suff);
addParameter(iP, 'SlopemIT', slopemITDefault, validfn_suff);
addParameter(iP, 'SlopehIT', slopehITDefault, validfn_suff);
addParameter(iP, 'GhbarIh', ghbarIhDefault, validfn_suff);
addParameter(iP, 'ShiftmIh', shiftmIhDefault, validfn_suff);
addParameter(iP, 'GkbarIKir', gkbarIKirDefault, validfn_suff);
addParameter(iP, 'GkbarIA', gkbarIADefault, validfn_suff);
addParameter(iP, 'GnabarINaP', gnabarINaPDefault, validfn_suff);

addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotInfFlag', plotInfFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotAllFlag', plotAllFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
v = iP.Results.v;
celsius = iP.Results.Celsius;
eK = iP.Results.Ek;
eNa = iP.Results.Ena;
eh = iP.Results.Eh;
caOut = iP.Results.CaOut;
caIn = iP.Results.CaIn;
pcabarIT = iP.Results.PcabarIT;
shiftmIT = iP.Results.ShiftmIT;
shifthIT = iP.Results.ShifthIT;
slopemIT = iP.Results.SlopemIT;
slopehIT = iP.Results.SlopehIT;
ghbarIh = iP.Results.GhbarIh;
shiftmIh = iP.Results.ShiftmIh;
gkbarIKir = iP.Results.GkbarIKir;
gkbarIA = iP.Results.GkbarIA;
gnabarINaP = iP.Results.GnabarINaP;
outFolder = iP.Results.OutFolder;
%suffices = iP.Results.Suffices;
plotInfFlag = iP.Results.PlotInfFlag;
plotAllFlag = iP.Results.PlotAllFlag;

% Construct voltage vector if not provided
if isempty(v)
    v = vMinDefault:vStepDefault:vMaxDefault;
else
    vMin = min(v);
    vMax = max(v);
end

% Compute IV curves for each current without plotting
[~, ~, ~, ~, ~, IMaxIT, IInitIT, IInfIT] = ...
    m3ha_compute_and_plot_IT2(v, 'Suffices', suffices, 'Celsius', celsius, ...
            'Cout', caOut, 'Cin', caIn, 'Pcabar', pcabarIT, ...
            'Shiftm', shiftmIT, 'Shifth', shifthIT, ...
            'Slopem', slopemIT, 'Slopeh', slopehIT, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);
[~, ~, ~, IMaxIh, IInfIh] = ...
    m3ha_compute_and_plot_Ih2(v, 'Suffices', suffices, 'Celsius', celsius, ...
            'Ghbar', ghbarIh, 'Erev', eh, 'Shiftm', shiftmIh, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);
[~, ~, IMaxIKir, IInfIKir] = ...
    m3ha_compute_and_plot_IKir2(v, 'Suffices', suffices, 'Celsius', celsius, ...
            'Gkbar', gkbarIKir, 'Ek', eK, ...
            'PlotInfFlag', false, 'PlotIVFlag', false);
[~, ~, ~, ~, ~, ~, ~, IMaxIA, IInitIA, IInfIA] = ...
    m3ha_compute_and_plot_IA2(v, 'Suffices', suffices, 'Celsius', celsius, ...
            'Gkbar', gkbarIA, 'Ek', eK, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);
[~, ~, ~, ~, IMaxINaP, IInitINaP, IInfINaP] = ...
    m3ha_compute_and_plot_INaP2(v, 'Suffices', suffices, 'Celsius', celsius, ...
            'GNabar', gnabarINaP, 'ENa', eNa, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);

% Compute total I-V curves
IAllInit = cellfun(@(a,b,c,d,f) a + b + c + d + f, IInitIT, IInfIh, IInfIKir, IInitIA, IInitINaP, 'UniformOutput', false);
IAllInf = cellfun(@(a,b,c,d,f) a + b + c + d + f, IInfIT, IInfIh, IInfIKir, IInfIA, IInfINaP, 'UniformOutput', false);

% Plot summary I-V curves
compute_and_plot_summary_IV(outFolder, suffices, v, IInitIT, IInfIT, IInfIh, ...
                IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP, IAllInit, IAllInf);

% Plot steady state I-V curves only
compute_and_plot_steady_IV(outFolder, suffices, v, IInfIT, IInfIh, ...
                IInfIKir, IInfIA, IInfINaP, IAllInf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hfig = compute_and_plot_summary_IV (outFolder, suffices, v, IInitIT, IInfIT, IInfIh, IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP, IAllInit, IAllInf)
linestyles = {'--', '-'};
hfig_all = figure('Visible', 'off');
xlim([min(v), max(v)]);
xlabel('Membrane potential (mV)');
ylabel('Current (mA/cm^2)');
title(['I-V relationships for ' strrep([suffices{:}], '_', '\_')]);
% set(gca, 'FontSize', 20);
for iSuffix = 1:numel(suffices)
	fprintf('Plotting summary I-V curve for %s ... \n\n', suffices{iSuffix}(2:end));
	hfig = figure('Visible', 'off');
	clf(hfig);
	hold on;
	cm = colormap(jet(5));
	plot(v, IInitIT{iSuffix}, '--', 'LineWidth', 2, 'Color', cm(1, :), ...
		'DisplayName', 'I_{T,init}');
	plot(v, IInfIT{iSuffix}, 'LineWidth', 2, 'Color', cm(1, :), ...
		'DisplayName', 'I_{T,\infty}');
	plot(v, IInfIh{iSuffix}, 'LineWidth', 2, 'Color', cm(2, :), ...
		'DisplayName', 'I_{h,\infty}');
	plot(v, IInfIKir{iSuffix}, 'LineWidth', 2, 'Color', cm(3, :), ...
		'DisplayName', 'I_{Kir,\infty}');
	plot(v, IInitIA{iSuffix}, '--', 'LineWidth', 2, 'Color', cm(4, :), ...
		'DisplayName', 'I_{A,init}');
	plot(v, IInfIA{iSuffix}, 'LineWidth', 2, 'Color', cm(4, :), ...
		'DisplayName', 'I_{A,\infty}');
	plot(v, IInitINaP{iSuffix}, '--', 'LineWidth', 2, 'Color', cm(5, :), ...
		'DisplayName', 'I_{NaP,init}');
	plot(v, IInfINaP{iSuffix}, 'LineWidth', 2, 'Color', cm(5, :), ...
		'DisplayName', 'I_{NaP,\infty}');
	plot(v, IAllInit{iSuffix}, '--', 'LineWidth', 2, 'Color', 'k', ...
		'DisplayName', 'I_{All,init}');
	plot(v, IAllInf{iSuffix}, 'LineWidth', 2, 'Color', 'k', ...
		'DisplayName', 'I_{All,\infty}');
	xlim([min(v), max(v)]);
	xlabel('Membrane potential (mV)');
	ylabel('Current (mA/cm^2)');
	title(['I-V relationships for ', strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	legend('location', 'northwest');
	% set(gca, 'FontSize', 20);
	saveas(hfig, fullfile(outFolder, ['I-V_all', suffices{iSuffix}]), 'png');

	set(0, 'CurrentFigure', hfig_all);
	hold on;
	cm2 = colormap(jet(8));
	plot(v, IInitIT{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(1, :), ...
		'DisplayName', ['I_{T,init}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfIT{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(2, :), ...
		'DisplayName', ['I_{T,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfIh{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(3, :), ...
		'DisplayName', ['I_{h,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfIKir{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(4, :), ...
		'DisplayName', ['I_{Kir,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInitIA{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(5, :), ...
		'DisplayName', ['I_{A,init}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfIA{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(6, :), ...
		'DisplayName', ['I_{A,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInitINaP{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(7, :), ...
		'DisplayName', ['I_{NaP,init}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfINaP{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm2(8, :), ...
		'DisplayName', ['I_{NaP,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IAllInit{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', 'k', ...
		'DisplayName', ['I_{All,init}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IAllInf{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', [0.1 0.1 0.1], ...
		'DisplayName', ['I_{All,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
end
legend('location', 'eastoutside');
saveas(hfig_all, fullfile(outFolder, ['I-V_all', suffices{:}]), 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hfig = compute_and_plot_steady_IV (outFolder, suffices, v, IInfIT, IInfIh, IInfIKir, IInfIA, IInfINaP, IAllInf)

linestyles = {'--', '-'};
hfig_all = figure('Visible', 'off');
xlim([min(v), max(v)]);
xlabel('Membrane potential (mV)');
ylabel('Current (mA/cm^2)');
title(['I-V relationships for ' strrep([suffices{:}], '_', '\_')]);
% set(gca, 'FontSize', 20);
for iSuffix = 1:numel(suffices)
	fprintf('Plotting steady state I-V curve for %s ... \n\n', suffices{iSuffix}(2:end));
	hfig = figure('Visible', 'off');
	clf(hfig);
	hold on;
	cm = colormap(jet(5));
	plot(v, IInfIT{iSuffix}, 'LineWidth', 2, 'Color', cm(1, :), ...
		'DisplayName', 'I_{T,\infty}');
	plot(v, IInfIh{iSuffix}, 'LineWidth', 2, 'Color', cm(2, :), ...
		'DisplayName', 'I_{h,\infty}');
	plot(v, IInfIKir{iSuffix}, 'LineWidth', 2, 'Color', cm(3, :), ...
		'DisplayName', 'I_{Kir,\infty}');
	plot(v, IInfIA{iSuffix}, 'LineWidth', 2, 'Color', cm(4, :), ...
		'DisplayName', 'I_{A,\infty}');
	plot(v, IInfINaP{iSuffix}, 'LineWidth', 2, 'Color', cm(5, :), ...
		'DisplayName', 'I_{NaP,\infty}');
	plot(v, IAllInf{iSuffix}, 'LineWidth', 2, 'Color', 'k', ...
		'DisplayName', 'I_{All,\infty}');
	xlim([min(v), max(v)]);
	xlabel('Membrane potential (mV)');
	ylabel('Current (mA/cm^2)');
	title(['I-V relationships for ', strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	legend('location', 'southeast');
	% set(gca, 'FontSize', 20);
	saveas(hfig, fullfile(outFolder, ['I-V_steady', suffices{iSuffix}]), 'png');

	set(0, 'CurrentFigure', hfig_all);
	hold on;
	plot(v, IInfIT{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm(1, :), ...
		'DisplayName', ['I_{T,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfIh{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm(2, :), ...
		'DisplayName', ['I_{h,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfIKir{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm(3, :), ...
		'DisplayName', ['I_{Kir,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfIA{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm(4, :), ...
		'DisplayName', ['I_{A,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IInfINaP{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', cm(5, :), ...
		'DisplayName', ['I_{NaP,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	plot(v, IAllInf{iSuffix}, 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'Color', 'k', ...
		'DisplayName', ['I_{All,\infty}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
end
legend('location', 'eastoutside');
saveas(hfig_all, fullfile(outFolder, ['I-V_steady', suffices{:}]), 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

IAllInit = IInitIT + IInfIh + IInfIKir + IInitIA + IInitINaP;
IAllInf =  IInfIT +  IInfIh + IInfIKir + IInfIA + IInfINaP;
hfig_all.Visible = 'off';
	hfig.Visible = 'off';
hfig_all.Visible = 'off';
	hfig.Visible = 'off';

%}

