function [v, minf, hinf, taum, tauh, IMax, IInit, IInf] = m3ha_compute_and_plot_IT2(varargin)
%% Plot activation/inactivation curves for the T-type calcium current
% Usage: [v, minf, hinf, taum, tauh, IMax, IInit, IInf] = m3ha_compute_and_plot_IT2(varargin)
%
% Arguments:    
%       v           (opt) voltage vector (x-values of plots)
%                   must be an increasing numeric vector
%                   default: -120:0.5:30 mV
%       varargin    - 'Celsius': temperature of simulation [degC]
%                   must be a numeric scalar
%                   default == 33 degC
%                   - 'Cout': extracellular [Ca++] [mM]
%                   must be a numeric scalar
%                   default == 2 mM
%                   - 'Cin': intracellular [Ca++] [mM]
%                   must be a numeric scalar
%                   default == 2.4e-4 mM
%                   - 'Pcabar': maximum Ca++ permeability [cm/s]
%                   must be a numeric scalar
%                   default == 0.2e-3 cm/s
%                   - 'Shiftm': depolarizing shift of activation curves [mV]
%                   must be a numeric scalar
%                   default == 1 mV
%                   - 'Shifth': depolarizing shift of inactivation curves [mV]
%                   must be a numeric scalar
%                   default == 1 mV
%                   - 'Slopem': scaling factor for slope of activation curves
%                   must be a numeric scalar
%                   default == 1
%                   - 'Slopeh': scaling factor for slope of inactivation curves
%                   must be a numeric scalar
%                   default == 1
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
%       cd/m3ha_compute_minf_IT.m
%       cd/m3ha_compute_hinf_IT.m
%       cd/m3ha_compute_taum_IT.m
%       cd/m3ha_compute_tauh_IT.m
%       /home/Matlab/Adams_Functions/compute_IMax_GHK.m
%
% File History:
% 2017-09-22 BT - Adapted from m3ha_compute_and_plot_IT.m
% 2017-10-04 BT - Inputs must be size == numel(suffixes), except celsius, cOut, and cIn
% 2017-10-20 BT - Plots suffixes against each other
% 2018-03-08 AL - Made figures invisible upon creation

%% Fixed values for IT
z = 2;                  % valence

%% Default voltage vector limits and steps
vMinDefault = -120;         % lowest voltage to plot [mV]
vMaxDefault = 30;           % highest voltage to plot [mV]
vStepDefault = 0.05;        % voltage step for plotting [mV]

%% Default parameters
celsiusDefault = 33;        % default temperature of simulation [degC]
cOutDefault = 2;            % default extracellular [Ca++] [mM]
cInDefault = 2.4e-4;        % default intracellular [Ca++] [mM]
pcabarDefault = 0.2e-3;     % default maximum Ca++ permeability [cm/s]
shiftmDefault = 1;  % default depolarizing shift of activation curves [mV]
shifthDefault = 1;  % default depolarizing shift of inactivation curves [mV]
slopemDefault = 1;  % default scaling factor for slope of activation curves
slopehDefault = 1;  % default scaling factor for slope of inactivation curves
outFolderDefault = pwd;     % default output directory
suffixesDefault = '';         % default suffixes for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves by default
plotTauFlagDefault = true;  % whether to plot time constant curves by default
plotIVFlagDefault = true;   % whether to plot I-V curves by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_IT2';

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
validfn_single = @(x) validateattributes(x, {'numeric'}, {'scalar'});	% celsius, cOut, cIn are single
%celsiusDefault = ones(1,numel(suffixes)) * celsiusDefault;
%cOutDefault = ones(1,numel(suffixes)) * cOutDefault;
%cInDefault = ones(1,numel(suffixes)) * cInDefault;
pcabarDefault = ones(1,numel(suffixes)) * pcabarDefault;
shiftmDefault = ones(1,numel(suffixes)) * shiftmDefault;
shifthDefault = ones(1,numel(suffixes)) * shifthDefault;
slopemDefault = ones(1,numel(suffixes)) * slopemDefault;
slopehDefault = ones(1,numel(suffixes)) * slopehDefault;
addParameter(iP, 'Celsius', celsiusDefault, validfn_single);
addParameter(iP, 'Cout', cOutDefault, validfn_single);
addParameter(iP, 'Cin', cInDefault, validfn_single);
addParameter(iP, 'Pcabar', pcabarDefault, validfn_suff);
addParameter(iP, 'Shiftm', shiftmDefault, validfn_suff);
addParameter(iP, 'Shifth', shifthDefault, validfn_suff);
addParameter(iP, 'Slopem', slopemDefault, validfn_suff);
addParameter(iP, 'Slopeh', slopehDefault, validfn_suff);

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
cOut = iP.Results.Cout;
cIn = iP.Results.Cin;
pcabar = iP.Results.Pcabar;
shiftm = iP.Results.Shiftm;
shifth = iP.Results.Shifth;
slopem = iP.Results.Slopem;
slopeh = iP.Results.Slopeh;
outFolder = iP.Results.OutFolder;
%suffixes = iP.Results.Suffixes;
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

% Compute the steady state value of the activation gating variable
for iSuffix = 1:numel(suffixes)
	minf{iSuffix} = m3ha_compute_minf_IT (v, shiftm(iSuffix), slopem(iSuffix));
end
minfsq = cellfun(@(x) x.^2, minf, 'UniformOutput', false);

% Compute the steady state value of the inactivation gating variable
for iSuffix = 1:numel(suffixes)
	hinf{iSuffix} = m3ha_compute_hinf_IT (v, shiftm(iSuffix), slopem(iSuffix));
end

% Compute the time constant for activation
for iSuffix = 1:numel(suffixes)
	taum{iSuffix} = m3ha_compute_taum_IT (v, shiftm(iSuffix), slopem(iSuffix), celsius);
end

% Compute the time constant for inactivation
for iSuffix = 1:numel(suffixes)
	tauh{iSuffix} = m3ha_compute_tauh_IT (v, shiftm(iSuffix), slopem(iSuffix), celsius);
end

% Compute maximum current possible
for iSuffix = 1:numel(suffixes)
	IMax{iSuffix} = compute_IMax_GHK (v, pcabar(iSuffix), z, celsius, cOut, cIn);
end
IInit = cellfun(@(x,y) x.*y, IMax, minfsq, 'UniformOutput', false);
IInf = cellfun(@(x,y,z) x.*y.*z, IMax, minfsq, hinf, 'UniformOutput', false);
%IInit = IMax .* minfsq;
%IInf = IMax .* minfsq .* hinf;

% Plot and save the steady state values
linestyles = {'--','-'}; % bef is --, aft is -
if plotInfFlag
	hfig1_all = figure('Visible', 'off');
	ylim([0, 1]);
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Steady state value');
	title(['Steady state values for activation/inactivation of IT for ' strrep([suffixes{:}], '_', '\_')]);
	% set(gca, 'FontSize', 15);
	for iSuffix = 1:numel(suffixes)
		hfig1 = figure('Visible', 'off');
		clf(hfig1);
		hold on;
		plot(v, minf{iSuffix}, 'b--', 'LineWidth', 2, 'DisplayName', 'm_{\infty}');
		plot(v, minfsq{iSuffix}, 'b', 'LineWidth', 2, 'DisplayName', 'm_{\infty}^2');
		plot(v, hinf{iSuffix}, 'r', 'LineWidth', 2, 'DisplayName', 'h_{\infty}');
		ylim([0, 1]);
		xlim([vMin, vMax]);
		xlabel('Membrane potential (mV)');
		ylabel('Steady state value');
		title(['Steady state values for activation/inactivation of IT for ' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		legend('location', 'northeast');
		% set(gca, 'FontSize', 15);
		saveas(hfig1, fullfile(outFolder, ['IT_minf_hinf', suffixes{iSuffix}]), 'png');

		set(0, 'CurrentFigure', hfig1_all);
		hold on;
		plot(v, minf{iSuffix}, 'b', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'DisplayName', ['m_{\infty}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		plot(v, minfsq{iSuffix}, 'c', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'DisplayName', ['m_{\infty}^2' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		plot(v, hinf{iSuffix}, 'r', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'DisplayName', ['h_{\infty}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
	end
	legend('location', 'eastoutside');
	saveas(hfig1_all, fullfile(outFolder, ['IT_minf_hinf' suffixes{:}]), 'png');
end

% Plot and save the time constants
if plotTauFlag
	hfig2_all = figure('Visible', 'off');
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Time constant (ms)');
	title(['Time constants for activation/inactivation of IT for ' strrep([suffixes{:}], '_', '\_')]);
	% set(gca, 'FontSize', 15);
	for iSuffix = 1:numel(suffixes)
		hfig2 = figure('Visible', 'off');
		clf(hfig2);
		hold on;
		plot(v, taum{iSuffix}, 'b', 'LineWidth', 2, 'DisplayName', '\tau_{m}');
		plot(v, tauh{iSuffix}, 'r', 'LineWidth', 2, 'DisplayName', '\tau_{h}');
		xlim([vMin, vMax]);
		xlabel('Membrane potential (mV)');
		ylabel('Time constant (ms)');
		title(['Time constants for activation/inactivation of IT for ' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		legend('location', 'northeast');
		% set(gca, 'FontSize', 15);
		saveas(hfig2, fullfile(outFolder, ['IT_taum_tauh', suffixes{iSuffix}]), 'png');

		set(0, 'CurrentFigure', hfig2_all);
		hold on;
		plot(v, taum{iSuffix}, 'b', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'DisplayName', ['\tau_{m}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		plot(v, tauh{iSuffix}, 'r', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, 'DisplayName', ['\tau_{h}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
	end
	legend('location', 'eastoutside');
	saveas(hfig2_all, fullfile(outFolder, ['IT_taum_tauh' suffixes{:}]), 'png');
end

% Plot and save the maximum current-voltage relationship
if plotIVFlag
	hfig3_all = figure('Visible', 'off');
	xlim([vMin, vMax]);
	xlabel('Membrane potential (mV)');
	ylabel('Current (mA/cm^2)');
	title(['I-V relationship of IT for ' strrep([suffixes{:}], '_', '\_')]);
	% set(gca, 'FontSize', 15);
	for iSuffix = 1:numel(suffixes)
		hfig3 = figure('Visible', 'off');
		clf(hfig3);
		hold on;
		plot(v, IMax{iSuffix}, 'k', 'LineWidth', 2, ...
		    'DisplayName', 'I_{max}');
		plot(v, IInit{iSuffix}, 'b--', 'LineWidth', 2, ...
		    'DisplayName', 'm_{\infty}^2I_{max}');
		plot(v, IInf{iSuffix}, 'b', 'LineWidth', 2, ...
		    'DisplayName', 'm_{\infty}^2h_{\infty}I_{max}');
		xlim([vMin, vMax]);
		xlabel('Membrane potential (mV)');
		ylabel('Current (mA/cm^2)');
		title(['I-V relationship of IT for ', strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		legend('location', 'southeast');
		% set(gca, 'FontSize', 15);
		saveas(hfig3, fullfile(outFolder, ['IT_I-V', suffixes{iSuffix}]), 'png');

		set(0, 'CurrentFigure', hfig3_all);
		hold on;
		plot(v, IMax{iSuffix}, 'k', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
			'DisplayName', ['I_{max}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		plot(v, IInit{iSuffix}, 'c', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
			'DisplayName', ['m_{\infty}^2I_{max}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
		plot(v, IInf{iSuffix}, 'b', 'LineStyle', linestyles{iSuffix}, 'LineWidth', 2, ...
			'DisplayName', ['m_{\infty}^2h_{\infty}I_{max}' strrep(suffixes{iSuffix}(2:end), '_', '\_')]);
	end
	legend('location', 'eastoutside');
	saveas(hfig3_all, fullfile(outFolder, ['IT_I-V' suffixes{:}]), 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

minf = m3ha_compute_minf_IT (v, shiftm, slopem);
hinf = m3ha_compute_hinf_IT (v, shiftm, slopem);
taum = m3ha_compute_taum_IT (v, shiftm, slopem, celsius);

% Compute the time constant for inactivation
tauh = m3ha_compute_tauh_IT (v, shiftm, slopem, celsius);

% Compute maximum current possible
IMax = compute_IMax_GHK (v, pcabar, z, celsius, cOut, cIn);
IInit = IMax .* minfsq;
IInf = IMax .* minfsq .* hinf;
red = [1 0 0];
for iSuffix = 1:numel(suffixes)
	red(iSuffix,:) = red(1,:) / numel(suffixes) * (numel(suffixes)-iSuffix+1);
end;
red = flipud(red);
blue = fliplr(red);
black = repmat(red(:,1), 1, 3);
black(end,:) = [0 0 0];

hfig1_all.Visible = 'off';
	hfig1.Visible = 'off';
hfig2_all.Visible = 'off';
	hfig2.Visible = 'off';
hfig3_all.Visible = 'off';
	hfig3.Visible = 'off';

%}

