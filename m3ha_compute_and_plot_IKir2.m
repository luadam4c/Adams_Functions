function [v, minf, IMax, IInf] = m3ha_compute_and_plot_IKir2(varargin)
%% Plot activation/inactivation curves for the T-type calcium current
% Usage: [v, minf, IMax, IInf] = m3ha_compute_and_plot_IKir2(varargin)
%
% Arguments:    
%       v           (opt) voltage vector (x-values of plots)
%                   must be an increasing numeric vector
%                   default: -120:0.5:30 mV
%       varargin    - 'Celsius': temperature of simulation [degC]
%                   must be a numeric scalar
%                   default == 33 degC
%                   - 'Gkbar': maximum conductance of IKir [S/cm2]
%                   must be a numeric scalar
%                   default == 2.0e-5 S/cm2
%                   - 'Ek': K+ reversal potential [mV]
%                   must be a numeric scalar
%                   default == -100 mV
%                   - 'OutFolder': output directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Suffices': suffices for figure names (includes underscore)
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'PlotInfFlag': whether to plot steady state curves
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotIVFlag': whether to plot I-V curves
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compare_and_plot_across_IC.m
%
% Requires:
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_minf_IKir.m
%
% File History:
% 2017-09-22 BT - Adapted from m3ha_compute_and_plot_IKir.m
% 2017-10-04 BT - Inputs must be size == numel(suffices), except celsius and eK
% 2018-03-08 AL - Made figures invisible upon creation

%% Voltage vector limits and steps
vMinDefault = -120;         % lowest voltage to plot [mV]
vMaxDefault = 30;           % highest voltage to plot [mV]
vStepDefault = 0.05;        % voltage step for plotting [mV]

%% Default parameters
celsiusDefault = 33;        % default temperature of simulation [degC]
gkbarDefault = 2.0e-5;      % default maximum conductance of IKir [S/cm2]
eKDefault = -100;           % default K+ reversal potential [mV]
outFolderDefault = pwd;     % default output directory
sufficesDefault = '';         % default suffices for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves by default
plotIVFlagDefault = true;   % whether to plot I-V curves by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_IKir2';

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
validfn_single = @(x) validateattributes(x, {'numeric'}, {'scalar'});	% celsius and Ek are single
%celsiusDefault = ones(1,numel(suffices)) * celsiusDefault;
gkbarDefault = ones(1,numel(suffices)) * gkbarDefault;
%eKDefault = ones(1,numel(suffices)) * eKDefault;
addParameter(iP, 'Celsius', celsiusDefault, validfn_single);
addParameter(iP, 'Gkbar', gkbarDefault, validfn_suff);
addParameter(iP, 'Ek', eKDefault, validfn_single);

addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotInfFlag', plotInfFlagDefault, ...
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
%suffices = iP.Results.Suffices;
plotInfFlag = iP.Results.PlotInfFlag;
plotIVFlag = iP.Results.PlotIVFlag;

% Construct voltage vector if not provided
if isempty(v)
    v = vMinDefault:vStepDefault:vMaxDefault;
else
    vMin = min(v);
    vMax = max(v);
end

% Compute the steady state values of the activation gating variable
minf = m3ha_compute_minf_IKir (v);

% Compute maximum current possible
for iSuffix = 1:numel(suffices)
	IMax{iSuffix} = gkbar(iSuffix) .* (v - eK);
end
IInf = cellfun(@(x) x .* minf, IMax, 'UniformOutput', false);

% Plot and save the steady state values
if plotInfFlag
	hfig1 = figure('Visible', 'off');
	clf(hfig1);
	hold on;
	plot(v, minf, 'b', 'LineWidth', 2, 'DisplayName', 'm_{\infty}');
	ylim([0, 1]);
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Steady state value');
	title('Steady state values for activation/inactivation of IKir');
	legend('location', 'northeast');
	% set(gca, 'FontSize', 15);
	saveas(hfig1, fullfile(outFolder, 'IKir_minf'), 'png');
end

% Plot and save the maximum current-voltage relationship
linestyles = {'--','-'}; % bef is --, aft is -
if plotIVFlag
	hfig3_all = figure('Visible', 'off');
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Current (mA/cm^2)');
	title(['I-V relationship of IKir for ' strrep([suffices{:}], '_', '\_')]);
	% set(gca, 'FontSize', 15);
	for iSuffix = 1:numel(suffices)
		hfig3 = figure('Visible', 'off');
		clf(hfig3);
		hold on;
		plot(v, IMax{iSuffix}, 'k', 'LineWidth', 2, ...
		    'DisplayName', 'I_{max}');
		plot(v, IInf{iSuffix}, 'b', 'LineWidth', 2, ...
		    'DisplayName', 'm_{\infty}I_{max}');
		xlim([vMin, vMax]);
		xlabel('Membrane potential (mV)');
		ylabel('Current (mA/cm^2)');
		title(['I-V relationship of IKir for ', strrep(suffices{iSuffix}(2:end), '_', '\_')]);
		legend('location', 'northwest');
		% set(gca, 'FontSize', 15);
		saveas(hfig3, fullfile(outFolder, ['IKir_I-V', suffices{iSuffix}]), 'png');

		set(0, 'CurrentFigure', hfig3_all);
		hold on;
		plot(v, IMax{iSuffix}, 'k', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
			'DisplayName', ['I_{max}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
		plot(v, IInf{iSuffix}, 'b', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
			'DisplayName', ['m_{\infty}I_{max}' strrep(suffices{iSuffix}(2:end), '_', '\_')]);
	end
	legend('location', 'eastoutside');
	saveas(hfig3_all, fullfile(outFolder, ['IKir_I-V' suffices{:}]), 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%IMax = gkbar .* (v - eK);
%IInf = IMax .* minf;

	hfig1.Visible = 'off';
	hfig3_all.Visible = 'off';
		hfig3.Visible = 'off';

%}
