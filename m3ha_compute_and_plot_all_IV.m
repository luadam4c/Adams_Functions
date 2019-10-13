function [v, IInitIT, IInfIT, IInfIh, IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP] = m3ha_compute_and_plot_all_IV(varargin)
%% Plot I-V curves of all currents together
% Usage: [v, IInitIT, IInfIT, IInfIh, IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP] = m3ha_compute_and_plot_all_IV(varargin)
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
%                   - 'Suffix': suffix for figure names (includes underscore)
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
% 2017-08-09 Adapted from m3ha_compare_and_plot_across_IC.m
% 2017-08-10 Added suffix to titles
% 2017-08-23 Added IAllInit adn IAllInf
% 2018-03-08 Made figures invisible

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
suffixDefault = '';         % default suffix for figure names
plotInfFlagDefault = true;  % whether to plot steady state curves only
plotAllFlagDefault = true;  % whether to plot all I-V curves by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compute_and_plot_all_IV';

% Add optional inputs to the Input Parser
addOptional(iP, 'v', [], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Celsius', celsiusDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Ek', eKDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Ena', eNaDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Eh', ehDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'CaOut', caOutDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'CaIn', caInDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PcabarIT', pcabarITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ShiftmIT', shiftmITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ShifthIT', shifthITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'SlopemIT', slopemITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'SlopehIT', slopehITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GhbarIh', ghbarIhDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ShiftmIh', shiftmIhDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GkbarIKir', gkbarIKirDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GkbarIA', gkbarIADefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GnabarINaP', gnabarINaPDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
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
suffix = iP.Results.Suffix;
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
    m3ha_compute_and_plot_IT(v, 'Suffix', suffix, 'Celsius', celsius, ...
            'Cout', caOut, 'Cin', caIn, 'Pcabar', pcabarIT, ...
            'Shiftm', shiftmIT, 'Shifth', shifthIT, ...
            'Slopem', slopemIT, 'Slopeh', slopehIT, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);
[~, ~, ~, IMaxIh, IInfIh] = ...
    m3ha_compute_and_plot_Ih(v, 'Suffix', suffix, 'Celsius', celsius, ...
            'Ghbar', ghbarIh, 'Erev', eh, 'Shiftm', shiftmIh, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);
[~, ~, IMaxIKir, IInfIKir] = ...
    m3ha_compute_and_plot_IKir(v, 'Suffix', suffix, 'Celsius', celsius, ...
            'Gkbar', gkbarIKir, 'Ek', eK, ...
            'PlotInfFlag', false, 'PlotIVFlag', false);
[~, ~, ~, ~, ~, ~, ~, IMaxIA, IInitIA, IInfIA] = ...
    m3ha_compute_and_plot_IA(v, 'Suffix', suffix, 'Celsius', celsius, ...
            'Gkbar', gkbarIA, 'Ek', eK, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);
[~, ~, ~, ~, IMaxINaP, IInitINaP, IInfINaP] = ...
    m3ha_compute_and_plot_INaP(v, 'Suffix', suffix, 'Celsius', celsius, ...
            'GNabar', gnabarINaP, 'ENa', eNa, ...
            'PlotInfFlag', false, 'PlotTauFlag', false, 'PlotIVFlag', false);

% Compute total I-V curves
IAllInit = IInitIT + IInfIh + IInfIKir + IInitIA + IInitINaP;
IAllInf = IInfIT + IInfIh + IInfIKir + IInfIA + IInfINaP;

% Plot summary I-V curves
compute_and_plot_summary_IV(outFolder, suffix, v, IInitIT, IInfIT, IInfIh, ...
                IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP, IAllInit, IAllInf);

% Plot steady state I-V curves only
compute_and_plot_steady_IV(outFolder, suffix, v, IInfIT, IInfIh, ...
                IInfIKir, IInfIA, IInfINaP, IAllInf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hfig = compute_and_plot_summary_IV (outFolder, suffix, v, IInitIT, IInfIT, IInfIh, IInfIKir, IInitIA, IInfIA, IInitINaP, IInfINaP, IAllInit, IAllInf)

fprintf('Plotting summary I-V curve for %s ... \n\n', suffix(2:end));
hfig = figure('Visible', 'off');
clf(hfig);
hold on;
cm = colormap((jet(5)));
plot(v, IInitIT, '--', 'LineWidth', 2, 'Color', cm(1, :), ...
    'DisplayName', 'I_{T,init}');
plot(v, IInfIT, 'LineWidth', 2, 'Color', cm(1, :), ...
    'DisplayName', 'I_{T,\infty}');
plot(v, IInfIh, 'LineWidth', 2, 'Color', cm(2, :), ...
    'DisplayName', 'I_{h,\infty}');
plot(v, IInfIKir, 'LineWidth', 2, 'Color', cm(3, :), ...
    'DisplayName', 'I_{Kir,\infty}');
plot(v, IInitIA, '--', 'LineWidth', 2, 'Color', cm(4, :), ...
    'DisplayName', 'I_{A,init}');
plot(v, IInfIA, 'LineWidth', 2, 'Color', cm(4, :), ...
    'DisplayName', 'I_{A,\infty}');
plot(v, IInitINaP, '--', 'LineWidth', 2, 'Color', cm(5, :), ...
    'DisplayName', 'I_{NaP,init}');
plot(v, IInfINaP, 'LineWidth', 2, 'Color', cm(5, :), ...
    'DisplayName', 'I_{NaP,\infty}');
plot(v, IAllInit, '--', 'LineWidth', 2, 'Color', 'k', ...
    'DisplayName', 'I_{All,init}');
plot(v, IAllInf, 'LineWidth', 2, 'Color', 'k', ...
    'DisplayName', 'I_{All,\infty}');
xlim([min(v), max(v)]);
xlabel('Membrane potential (mV)');
ylabel('Current (mA/cm^2)');
title(['I-V relationships for ', strrep(suffix(2:end), '_', '\_')]);
legend('location', 'northwest');
set(gca, 'FontSize', 20);
saveas(hfig, fullfile(outFolder, ['I-V_all', suffix]), 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hfig = compute_and_plot_steady_IV (outFolder, suffix, v, IInfIT, IInfIh, IInfIKir, IInfIA, IInfINaP, IAllInf)

fprintf('Plotting steady state I-V curve for %s ... \n\n', suffix(2:end));
hfig = figure('Visible', 'off');
clf(hfig);
hold on;
cm = colormap((jet(5)));
plot(v, IInfIT, 'LineWidth', 2, 'Color', cm(1, :), ...
    'DisplayName', 'I_{T,\infty}');
plot(v, IInfIh, 'LineWidth', 2, 'Color', cm(2, :), ...
    'DisplayName', 'I_{h,\infty}');
plot(v, IInfIKir, 'LineWidth', 2, 'Color', cm(3, :), ...
    'DisplayName', 'I_{Kir,\infty}');
plot(v, IInfIA, 'LineWidth', 2, 'Color', cm(4, :), ...
    'DisplayName', 'I_{A,\infty}');
plot(v, IInfINaP, 'LineWidth', 2, 'Color', cm(5, :), ...
    'DisplayName', 'I_{NaP,\infty}');
plot(v, IAllInf, 'LineWidth', 2, 'Color', 'k', ...
    'DisplayName', 'I_{All,\infty}');
xlim([min(v), max(v)]);
xlabel('Membrane potential (mV)');
ylabel('Current (mA/cm^2)');
title(['I-V relationships for ', strrep(suffix(2:end), '_', '\_')]);
legend('location', 'southeast');
set(gca, 'FontSize', 20);
saveas(hfig, fullfile(outFolder, ['I-V_steady', suffix]), 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

