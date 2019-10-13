function [v, minf, IMax, IInf] = m3ha_compute_and_plot_IKir(varargin)
%% Plot activation/inactivation curves for the T-type calcium current
% Usage: [v, minf, IMax, IInf] = m3ha_compute_and_plot_IKir(varargin)
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
%                   - 'Suffix': suffix for figure names (includes underscore)
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
% 2017-08-06 Adapted from plot_IA.m
% 2017-08-09 Added v as an optional argument
% 2017-08-10 Added suffix to titles
% 2018-03-08 Made figures invisible

%% Voltage vector limits and steps
vMinDefault = -120;         % lowest voltage to plot [mV]
vMaxDefault = 30;           % highest voltage to plot [mV]
vStepDefault = 0.05;        % voltage step for plotting [mV]

%% Default parameters
celsiusDefault = 33;        % default temperature of simulation [degC]
gkbarDefault = 2.0e-5;      % default maximum conductance of IKir [S/cm2]
eKDefault = -100;           % default K+ reversal potential [mV]
outFolderDefault = pwd;     % default output directory
suffixDefault = '';         % default suffix for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves by default
plotIVFlagDefault = true;   % whether to plot I-V curves by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_IKir';

% Add optional inputs to the Input Parser
addOptional(iP, 'v', [], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Celsius', celsiusDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Gkbar', gkbarDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Ek', eKDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
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
suffix = iP.Results.Suffix;
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
IMax = gkbar .* (v - eK);
IInf = IMax .* minf;

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
    set(gca, 'FontSize', 15);
    saveas(hfig1, fullfile(outFolder, ['IKir_minf', suffix]), 'png');
end

% Plot and save the maximum current-voltage relationship
if plotIVFlag
    hfig3 = figure('Visible', 'off');
    clf(hfig3);
    hold on;
    plot(v, IMax, 'k', 'LineWidth', 2, ...
        'DisplayName', 'I_{max}');
    plot(v, IInf, 'b', 'LineWidth', 2, ...
        'DisplayName', 'm_{\infty}I_{max}');
    xlim([vMin, vMax]);
    xlabel('Membrane potential (mV)');
    ylabel('Current (mA/cm^2)');
    title(['I-V relationship of IKir for ', strrep(suffix(2:end), '_', '\_')]);
    legend('location', 'northwest');
    set(gca, 'FontSize', 15);
    saveas(hfig3, fullfile(outFolder, ['IKir_I-V', suffix]), 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

