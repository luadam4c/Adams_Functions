function [v, m1inf, m2inf, hinf, taum, tauh1, tauh2, IMax, IInit, IInf] = m3ha_compute_and_plot_IA2(varargin)
%% Plot activation/inactivation curves for the T-type calcium current
% Usage: [v, m1inf, m2inf, hinf, taum, tauh1, tauh2, IMax, IInit, IInf] = m3ha_compute_and_plot_IA2(varargin)
%
% Arguments:    
%       v           (opt) voltage vector (x-values of plots)
%                   must be an increasing numeric vector
%                   default: -120:0.5:30 mV
%       varargin    - 'Celsius': temperature of simulation [degC]
%                   must be a numeric scalar
%                   default == 33 degC
%                   - 'Gkbar': maximum conductance of IA [S/cm2]
%                   must be a numeric scalar
%                   default == 5.5e-3 S/cm2
%                   - 'Ek': K+ reversal potential [mV]
%                   must be a numeric scalar
%                   default == -100 mV
%                   - 'OutFolder': output directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Suffixes': suffixes for figure names (includes underscore)
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'PlotInfFlag': whether to plot steady state curves
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotTauFlag': whether to plot time constant curves
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotIVFlag': whether to plot I-V curves
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Used by:    
%       cd/m3ha_compare_and_plot_across_IC.m
%
% Requires:
%       cd/m3ha_compute_m1inf_IA.m
%       cd/m3ha_compute_m2inf_IA.m
%       cd/m3ha_compute_hinf_IA.m
%       cd/m3ha_compute_taum_IA.m
%       cd/m3ha_compute_tauh1_IA.m
%       cd/m3ha_compute_tauh2_IA.m
%
% File History:
% 2017-09-22 BT - Adapted from m3ha_compute_and_plot_IA.m
% 2017-10-04 BT - Inputs must be size == numel(suffixes), except celsius and eK
% 2018-03-08 AL - Made figures invisible upon creation

%% Voltage vector limits and steps
vMinDefault = -120;         % lowest voltage to plot [mV]
vMaxDefault = 30;           % highest voltage to plot [mV]
vStepDefault = 0.05;        % voltage step for plotting [mV]

%% Default parameters
celsiusDefault = 33;        % default temperature of simulation [degC]
gkbarDefault = 5.5e-3;      % default maximum conductance of IA [S/cm2]
eKDefault = -100;           % default K+ reversal potential [mV]
outFolderDefault = pwd;     % default output directory
suffixesDefault = '';         % default suffix for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves by default
plotTauFlagDefault = true;  % whether to plot time constant curves by default
plotIVFlagDefault = true;   % whether to plot I-V curves by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_IA2';

% Check if Suffixes is scalar or cell array
addParameter(iP, 'Suffixes', suffixesDefault, ...
	@(x) validateattributes(x, {'char', 'string', 'cell'}, {'nonempty'}));
potSuffixesInd = find(strcmp(varargin, 'Suffixes'));	% potential index for Suffixes string
if ~isempty(potSuffixesInd)
	parse(iP, varargin{potSuffixesInd}, varargin{potSuffixesInd+1});	% parse suffixes if exists
else
	parse(iP);
end
suffixes = iP.Results.Suffixes;
if ~iscell(suffixes) & numel(suffixes) == 1	% if scalar text, convert to cell
	suffixes = {suffixes};
end

% Add optional inputs to the Input Parser
addOptional(iP, 'v', [], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}));

% Add parameter-value pairs to the Input Parser
validfn_suff = @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', numel(suffixes)});	% size match to numel(suffixes)
validfn_single = @(x) validateattributes(x, {'numeric'}, {'scalar'});	% celsius and Ek are single
%celsiusDefault = ones(1,numel(suffixes)) * celsiusDefault;	% change default values to size of suffixes
gkbarDefault = ones(1,numel(suffixes)) * gkbarDefault;
%eKDefault = ones(1,numel(suffixes)) * eKDefault;
addParameter(iP, 'Celsius', celsiusDefault, validfn_single);
addParameter(iP, 'Gkbar', gkbarDefault, validfn_suff);
addParameter(iP, 'Ek', eKDefault, validfn_single);

addParameter(iP, 'OutFolder', outFolderDefault, ...
	@(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotInfFlag', plotInfFlagDefault, ...
	@(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotTauFlag', plotTauFlagDefault, ...
	@(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotIVFlag', plotIVFlagDefault, ...
	@(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
v = iP.Results.v;
celsius = iP.Results.Celsius;
gkbar = iP.Results.Gkbar;
eK = iP.Results.Ek;
outFolder = iP.Results.OutFolder;
%suffix = iP.Results.Suffix;
plotInfFlag = iP.Results.PlotInfFlag;
plotTauFlag = iP.Results.PlotTauFlag;
plotIVFlag = iP.Results.PlotIVFlag;

% Construct voltage vector if not provided
if isempty(v)
    v = vMinDefault:vStepDefault:vMaxDefault;
else
    vMin = min(v);
    vMax = max(v);
end

% Compute the steady state values of the activation gating variables
m1inf = m3ha_compute_m1inf_IA (v);
m1infquad = m1inf.^4;
m2inf = m3ha_compute_m2inf_IA (v);
m2infquad = m2inf.^4;

% Compute the steady state value of the inactivation gating variables
hinf = m3ha_compute_hinf_IA (v);

% Compute the time constant for activation
taum = m3ha_compute_taum_IA (v, celsius);

% Compute the time constants for inactivation
tauh1 = m3ha_compute_tauh1_IA (v, celsius);
tauh2 = m3ha_compute_tauh2_IA (v, celsius);

% Compute maximum current possible
for iSuffix = 1:numel(suffixes)
	IMax{iSuffix} = gkbar(iSuffix) .* (v - eK);
	IInit{iSuffix} = IMax{iSuffix} .* (0.6 * m1infquad + 0.4 * m2infquad);
	IInf{iSuffix} = IMax{iSuffix} .* (0.6 * m1infquad .* hinf + 0.4 * m2infquad .* hinf);
end

% Plot and save the steady state values
if plotInfFlag
	hfig1 = figure('Visible', 'off');
	clf(hfig1);
	hold on;
	plot(v, m1inf, 'b--', 'LineWidth', 2, 'DisplayName', 'm_{1,\infty}');
	plot(v, m1infquad, 'b', 'LineWidth', 2, 'DisplayName', 'm_{1,\infty}^2');
	plot(v, m2inf, 'c--', 'LineWidth', 2, 'DisplayName', 'm_{2,\infty}');
	plot(v, m2infquad, 'c', 'LineWidth', 2, 'DisplayName', 'm_{2,\infty}^2');
	plot(v, hinf, 'r', 'LineWidth', 2, 'DisplayName', 'h_{\infty}');
	ylim([0, 1]);
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Steady state value');
	title('Steady state values for activation/inactivation of IA');
	legend('location', 'northeast');
	% set(gca, 'FontSize', 15);
	saveas(hfig1, fullfile(outFolder, 'IA_m1inf_m2inf_hinf'), 'png');
end

% Plot and save the time constants
linestyles = {'--','-'}; % bef is --, aft is -
if plotTauFlag
	hfig2 = figure('Visible', 'off');
	clf(hfig2);
	hold on;
	plot(v, taum, 'b', 'LineWidth', 2, 'DisplayName', '\tau_{m}');
	plot(v, tauh1, 'r', 'LineWidth', 2, 'DisplayName', '\tau_{h_1}');
	plot(v, tauh2, 'm', 'LineWidth', 2, 'DisplayName', '\tau_{h_2}');
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Time constant (ms)');
	title('Time constants for activation/inactivation of IA');
	legend('location', 'northeast');
	% set(gca, 'FontSize', 15);
	saveas(hfig2, fullfile(outFolder, 'IA_taum_tauh1_tauh2'), 'png');
end

if plotIVFlag
	hfig3_all = figure('Visible', 'off');
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Current (mA/cm^2)');
	title(['I-V relationship of IA for ' strrep([suffixes{:}], '_', '\_')]);
	% set(gca, 'FontSize', 15);
	for iSuffix = 1:numel(suffixes)
		% Plot and save the maximum current-voltage relationship
		hfig3 = figure('Visible', 'off');
		clf(hfig3);
		hold on;
		plot(v, IMax{iSuffix}, 'k', 'LineWidth', 2, ...
		    'DisplayName', 'I_{max}');
		plot(v, IInit{iSuffix}, 'b--', 'LineWidth', 2, ...
		    'DisplayName', '(0.6m_{1,\infty}^4+0.4m_{1,\infty}^4)I_{max}');
		plot(v, IInf{iSuffix}, 'b', 'LineWidth', 2, ...
		    'DisplayName', '(0.6m_{1,\infty}^4+0.4m_{1,\infty}^4)h_{\infty}I_{max}');
		xlim([vMin, vMax]);
		xlabel('Membrane potential (mV)');
		ylabel('Current (mA/cm^2)');
		title(['I-V relationship of IA for ', strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		legend('location', 'northwest');
		% set(gca, 'FontSize', 15);
		saveas(hfig3, fullfile(outFolder, ['IA_I-V', suffixes{iSuffix}]), 'png');
		
		set(0, 'CurrentFigure', hfig3_all);
		hold on;
		plot(v, IMax{iSuffix}, 'k', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
		    'DisplayName', ['I_{max}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		plot(v, IInit{iSuffix}, 'b', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
		    'DisplayName', ['(0.6m_{1,\infty}^4+0.4m_{1,\infty}^4)I_{max}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		plot(v, IInf{iSuffix}, 'c', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
		    'DisplayName', ['(0.6m_{1,\infty}^4+0.4m_{1,\infty}^4)h_{\infty}I_{max}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
	end
	legend('location', 'eastoutside');
	saveas(hfig3_all, fullfile(outFolder, ['IA_I-V' suffixes{:}]), 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

    @(x) assert((iscell(x) && (min(cellfun(@ischar, x)) ...
                                || min(cellfun(@isstring, x)))) || ...
				(isscalar(x) && ischar(x)), ...
        'suffixes must be a cell array of strings/character arrays or scalar text!'));
% IMax = gkbar .* (v - eK);
% IInit = IMax .* (0.6 * m1infquad + 0.4 * m2infquad);
% IInf = IMax .* (0.6 * m1infquad .* hinf + 0.4 * m2infquad .* hinf);

	hfig1.Visible = 'off';
	hfig2.Visible = 'off';
	hfig3_all.Visible = 'off';
		hfig3.Visible = 'off';

%}

