function [v, minf, taum, IMax, IInf] = m3ha_compute_and_plot_Ih(varargin)
%% Plot activation/inactivation curves for the T-type calcium current
% Usage: [v, minf, taum, IMax, IInf] = m3ha_compute_and_plot_Ih(varargin)
%
% Arguments:    
%       v           (opt) voltage vector (x-values of plots)
%                   must be an increasing numeric vector
%                   default: -120:0.5:30 mV
%       varargin    - 'Celsius': temperature of simulation [degC]
%                   must be a numeric scalar
%                   default == 33 degC
%                   - 'Ghbar': maximum conductance of Ih [S/cm2]
%                   must be a numeric scalar
%                   default == 2.2e-5 S/cm2
%                   - 'Erev': reversal potential of Ih [mV]
%                   must be a numeric scalar
%                   default == -43 mV
%                   - 'Shiftm': depolarizing shift of activation curves [mV]
%                   must be a numeric scalar
%                   default == 0 mV
%                   - 'OutFolder': output directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Suffix': suffix for figure names (includes underscore)
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
%       cd/m3ha_compute_minf_Ih.m
%       cd/m3ha_compute_taum_Ih.m
%
% File History:
% 2017-08-06 Created
% 2017-08-09 Added v as an optional argument
% 2017-08-10 Added suffix to titles
% 2018-03-08 Made figures invisible


%% Voltage vector limits and steps
vMinDefault = -120;            % lowest voltage to plot [mV]
vMaxDefault = 30;              % highest voltage to plot [mV]
vStepDefault = 0.05;           % voltage step for plotting [mV]

%% Default parameters
celsiusDefault = 33;        % default temperature of simulation [degC]
ghbarDefault = 2.2e-5;      % default maximum conductance of Ih [S/cm2]
eRevDefault = -43;          % default reversal potential of Ih [mV]
shiftmDefault = 0;          % default depolarizing shift of activation curves [mV]
outFolderDefault = pwd;     % default output directory
suffixDefault = '';         % default suffix for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves by default
plotTauFlagDefault = true;  % whether to plot time constant curves by default
plotIVFlagDefault = true;   % whether to plot I-V curves by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_Ih';

% Add optional inputs to the Input Parser
addOptional(iP, 'v', [], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Celsius', celsiusDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Ghbar', ghbarDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Erev', eRevDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Shiftm', shiftmDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
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
ghbar = iP.Results.Ghbar;
eRev = iP.Results.Erev;
shiftm = iP.Results.Shiftm;
outFolder = iP.Results.OutFolder;
suffix = iP.Results.Suffix;
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
minf = m3ha_compute_minf_Ih (v, shiftm);

% Compute the time constant for activation
taum = m3ha_compute_taum_Ih (v, shiftm, celsius);

% Compute maximum current possible
IMax = ghbar .* (v - eRev);
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
    title('Steady state values for activation of Ih');
    legend('location', 'northeast');
    set(gca, 'FontSize', 15);
    saveas(hfig1, fullfile(outFolder, ['Ih_minf', suffix]), 'png');
end

% Plot and save the time constants
if plotTauFlag
    hfig2 = figure('Visible', 'off');
    clf(hfig2);
    hold on;
    plot(v, taum, 'b', 'LineWidth', 2, 'DisplayName', '\tau_{m}');
    xlim([vMin, vMax]);
    xlabel('Membrane potential (mV)');
    ylabel('Time constant (ms)');
    title('Time constants for activation of Ih');
    legend('location', 'northeast');
    set(gca, 'FontSize', 15);
    saveas(hfig2, fullfile(outFolder, ['Ih_taum', suffix]), 'png');
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
    title(['I-V relationship of Ih for ', strrep(suffix(2:end), '_', '\_')]);
    legend('location', 'southeast');
    set(gca, 'FontSize', 15);
    saveas(hfig3, fullfile(outFolder, ['Ih_I-V', suffix]), 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

