function [v, minf, hinf, taum, tauh, IMax, IInit, IInf] = ...
                m3ha_compute_and_plot_IT (varargin)
%% Plot activation/inactivation curves for the T-type calcium current
% Usage: [v, minf, hinf, taum, tauh, IMax, IInit, IInf] = ...
%               m3ha_compute_and_plot_IT (varargin)
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
%       cd/m3ha_compute_minf_IT.m
%       cd/m3ha_compute_hinf_IT.m
%       cd/m3ha_compute_taum_IT.m
%       cd/m3ha_compute_tauh_IT.m
%       cd/update_figure_for_corel.m
%       /home/Matlab/Adams_Functions/compute_IMax_GHK.m
%
% File History:
% 2017-08-05 Created
% 2017-08-06 Added compute_IMax_GHK
% 2017-08-09 Added v as an optional argument
% 2017-08-10 Added suffix to titles
% 2018-03-08 Made figures invisible
% 2020-07-09 Now uses update_figure_for_corel.m

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
suffixDefault = '';         % default suffix for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves by default
plotTauFlagDefault = true;  % whether to plot time constant curves by default
plotIVFlagDefault = true;   % whether to plot I-V curves by default

% TODO:
figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_IT';

% Add optional inputs to the Input Parser
addOptional(iP, 'v', [], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Celsius', celsiusDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Cout', cOutDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Cin', cInDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Pcabar', pcabarDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Shiftm', shiftmDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Shifth', shifthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Slopem', slopemDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Slopeh', slopehDefault, ...
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
cOut = iP.Results.Cout;
cIn = iP.Results.Cin;
pcabar = iP.Results.Pcabar;
shiftm = iP.Results.Shiftm;
shifth = iP.Results.Shifth;
slopem = iP.Results.Slopem;
slopeh = iP.Results.Slopeh;
outFolder = iP.Results.OutFolder;
suffix = iP.Results.Suffix;
plotInfFlag = iP.Results.PlotInfFlag;
plotTauFlag = iP.Results.PlotTauFlag;
plotIVFlag = iP.Results.PlotIVFlag;

% Construct voltage vector if not provided
if isempty(v)
    v = vMinDefault:vStepDefault:vMaxDefault;
end

% Compute the minimum and maximum voltage
vMin = min(v);
vMax = max(v);

% Compute the steady state value of the activation gating variable
minf = m3ha_compute_minf_IT (v, shiftm, slopem);
minfsq = minf.^2;

% Compute the steady state value of the inactivation gating variable
hinf = m3ha_compute_hinf_IT (v, shiftm, slopem);

% Compute the time constant for activation
taum = m3ha_compute_taum_IT (v, shiftm, slopem, celsius);

% Compute the time constant for inactivation
tauh = m3ha_compute_tauh_IT (v, shiftm, slopem, celsius);

% Compute maximum current possible
IMax = compute_IMax_GHK (v, pcabar, z, celsius, cOut, cIn);
IInit = IMax .* minfsq;
IInf = IMax .* minfsq .* hinf;

% Compute maximum open probability at steady state
openProbability = minf .^ 2 .* hinf;
maxOpenProbability = max(openProbability);
fprintf("maxOpenProbability = %g\n", maxOpenProbability);

% Plot and save the steady state values
if plotInfFlag
    hfig1 = set_figure_properties('AlwaysNew', true);
    clf(hfig1);
    hold on;
    p1 = plot(v, minf, 'b--', 'LineWidth', 2, 'DisplayName', 'm_{\infty}');
    p2 = plot(v, minfsq, 'b', 'LineWidth', 2, 'DisplayName', 'm_{\infty}^2');
    p3 = plot(v, hinf, 'r', 'LineWidth', 2, 'DisplayName', 'h_{\infty}');
    ylim([0, 1]);
    xlim([vMin, vMax]);
    xlabel('Membrane potential (mV)');
    ylabel('Steady state value');
    % title('Steady state values for activation/inactivation of IT');
    legend('location', 'southwest');
    update_figure_for_corel(hfig1, 'LabelsFontSize', 15, ...
                            'AxisFontSize', 15, 'TextFontSize', 15);
    save_all_figtypes(hfig1, ...
        fullfile(outFolder, ['IT_minf_hinf', suffix]), figTypes);

    delete(p1);
    update_figure_for_corel(hfig1, 'RemoveLegends', true, ...
                            'LabelsFontSize', 15, ...
                            'AxisFontSize', 15, 'TextFontSize', 15);
    save_all_figtypes(hfig1, ...
        fullfile(outFolder, ['IT_minf_hinf_temp', suffix]), figTypes);
end

% Plot and save the time constants
if plotTauFlag
    hfig2 = set_figure_properties('AlwaysNew', true);
    clf(hfig2);
    hold on;
    plot(v, taum, 'b', 'LineWidth', 2, 'DisplayName', '\tau_{m}');
    plot(v, tauh, 'r', 'LineWidth', 2, 'DisplayName', '\tau_{h}');
    xlim([vMin, vMax]);
    xlabel('Membrane potential (mV)');
    ylabel('Time constant (ms)');
    % title('Time constants for activation/inactivation of IT');
    legend('location', 'northeast');
    set(gca, 'FontSize', 15);
    save_all_figtypes(hfig2, ...
        fullfile(outFolder, ['IT_taum_tauh', suffix]), figTypes);
end

% Plot and save the maximum current-voltage relationship
if plotIVFlag
    hfig3 = set_figure_properties('AlwaysNew', true);
    clf(hfig3);
    hold on;
    plot(v, IMax, 'k', 'LineWidth', 2, ...
        'DisplayName', 'I_{max}');
    plot(v, IInit, 'b--', 'LineWidth', 2, ...
        'DisplayName', 'm_{\infty}^2I_{max}');
    plot(v, IInf, 'b', 'LineWidth', 2, ...
        'DisplayName', 'm_{\infty}^2h_{\infty}I_{max}');
    xlim([vMin, vMax]);
    xlabel('Membrane potential (mV)');
    ylabel('Current (mA/cm^2)');
    % title(['I-V relationship of IT for ', strrep(suffix(2:end), '_', '\_')]);
    legend('location', 'southeast');
    set(gca, 'FontSize', 15);
    save_all_figtypes(hfig3, ...
        fullfile(outFolder, ['IT_I-V', suffix]), figTypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

