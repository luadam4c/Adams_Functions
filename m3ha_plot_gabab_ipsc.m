function varargout = m3ha_plot_gabab_ipsc (outFolder, varargin)
%% m3ha_plot_gabab_ipsc.m
% Requires:
%       cd/compute_gabab_conductance.m
%       cd/logscale.m
% TODO
%       cd/array_fun.m
%       cd/compute_time_constant.m
%       cd/create_labels_from_numbers.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/plot_traces.m
%
% Used by:
%       cd/m3ha_plot_figure05.m
%

% 2020-01-03 Moved from plot_gababipsc.m

%% Hard-coded parameters
siMs = 0.1;
tStart = 0;
tStop = 8000;
ipscStart = 0;                                              % (ms)
ampOrig = [32.00; 48.00; 17.76; 12.64];                     % (nS)
tauRiseOrig = [52.00; 52.00; 38.63; 39.88];                 % (ms)
tauFallFastOrig = [90.10; 90.10; 273.40; 65.80];            % (ms)
tauFallSlowOrig = [1073.20; 1073.20; 1022.00; 2600.00];     % (ms)
weight = [0.952; 0.952; 0.775; 0.629];

figTypes = {'png', 'epsc2'};
colorMap = [];

plotOriginal = false;
plotVaryWeight = false;
plotVaryTau = true;
plotVaryTauFall = false; %true;
plotVaryShapeOld = false;
plotVaryDualtoGAT3 = false;
plotVaryDualtoGAT3toGAT1 = false; %true;
plotOldFigures = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Create a time vector in ms
tVec = transpose(tStart:siMs:tStop);

%% Extract values for the Control condition
ampCon = ampOrig(1);
tauRiseCon = tauRiseOrig(1);
tauFallFastCon = tauFallFastOrig(1);
tauFallSlowCon = tauFallSlowOrig(1);
weightCon = weight(1);

%% Extract values for the GAT1 blockade condition
ampGat1 = ampOrig(2);
tauRiseGat1 = tauRiseOrig(2);
tauFallFastGat1 = tauFallFastOrig(2);
tauFallSlowGat1 = tauFallSlowOrig(2);
weightGat1 = weight(2);

%% Extract values for the GAT3 blockade condition
ampGat3 = ampOrig(3);
tauRiseGat3 = tauRiseOrig(3);
tauFallFastGat3 = tauFallFastOrig(3);
tauFallSlowGat3 = tauFallSlowOrig(3);
weightGat3 = weight(3);

%% Extract values for the Dual blockade condition
ampDual = ampOrig(4);
tauRiseDual = tauRiseOrig(4);
tauFallFastDual = tauFallFastOrig(4);
tauFallSlowDual = tauFallSlowOrig(4);
weightDual = weight(4);

%% Compute area under the curves
aucCon = compute_ipsc_area(tVec, ipscStart, ampCon, tauRiseCon, ...
                            tauFallFastCon, tauFallSlowCon, weightCon);
fprintf('aucCon = %g\n', aucCon);
aucGat1 = compute_ipsc_area(tVec, ipscStart, ampGat1, tauRiseGat1, ...
                            tauFallFastGat1, tauFallSlowGat1, weightGat1);
fprintf('aucGat1 = %g\n', aucGat1);
aucGat3 = compute_ipsc_area(tVec, ipscStart, ampGat3, tauRiseGat3, ...
                            tauFallFastGat3, tauFallSlowGat3, weightGat3);
% TODO for SHINSHIN: Create print_value.m
fprintf('aucGat3 = %g\n', aucGat3);
aucDual = compute_ipsc_area(tVec, ipscStart, ampDual, tauRiseDual, ...
                            tauFallFastDual, tauFallSlowDual, weightDual);
fprintf('aucDual = %g\n', aucDual);

%% Compute the GAT1 or GAT3 amplitude necessary to reach the same area 
%   under the curve as Dual blockade
%   and compute area under the curves for verification
ampGat1New = ampGat1 .* aucDual / aucGat1;
aucGat1New = compute_ipsc_area(tVec, ipscStart, ampGat1New, tauRiseGat1, ...
                            tauFallFastGat1, tauFallSlowGat1, weightGat1);
fprintf('aucGat1New = %g\n', aucGat1New);

ampGat3New = ampGat3 .* aucDual / aucGat3;
aucGat3New = compute_ipsc_area(tVec, ipscStart, ampGat3New, tauRiseGat3, ...
                            tauFallFastGat3, tauFallSlowGat3, weightGat3);
fprintf('aucGat3New = %g\n', aucGat3New);

%% Plot original GABAB IPSC conductances from Christine's thesis & old network model
if plotOriginal
    fig = set_figure_properties('AlwaysNew', true);
    figPathBase = fullfile(outFolder, 'gababipsc_original');
    defaultPosition = get(fig, 'Position');
    pos(1,:) = defaultPosition;
    pos(1, 1) = defaultPosition(1, 1) - defaultPosition(1, 3);
    set(fig, 'Position', pos(1,:));
    gVecs = compute_gabab_conductance(tVec, ipscStart, ampOrig, tauRiseOrig, ...
                                    tauFallFastOrig, tauFallSlowOrig, weight);
    if min(min(max(gVecs, 1))) ~= 1; error('Code error!'); end
    tauEmpirical = siMs .* compute_time_constant(gVecs, 'DecayMethod', 'empirical');
    plot_conductance(tVec, gVecs, colorMap);
    % legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
    legend(create_labels_from_numbers(tauEmpirical, 'Prefix', 'tau = '));
    title('GABAB IPSC conductances from thesis & old network model')
    save_all_figtypes(fig, figPathBase, figTypes);
end

%% Plot as weight is varied
if plotVaryWeight
    fig = set_figure_properties('AlwaysNew', true);
    figPathBase = fullfile(outFolder, 'gababipsc_dual_vary_weight');
    weightTest = 0:0.1:1;
    gVecs = compute_gabab_conductance(tVec, ipscStart, ampDual, tauRiseDual, ...
                                tauFallFastDual, tauFallSlowDual, weightTest);
    plot_conductance(tVec, gVecs, colorMap);
    legend(create_labels_from_numbers(weightTest, 'Prefix', 'weight = '));
    save_all_figtypes(fig, figPathBase, figTypes);
end

%% Plot as all taus are varied with the same
%   relative value and weights
if plotVaryTau
    fig = set_figure_properties('AlwaysNew', true);
    figPathBase = fullfile(outFolder, 'gababipsc_dual_vary_tau');
    scaleFactors = 0.1:0.1:1;
    tauRiseTest = tauRiseDual * scaleFactors;
    tauFallFastTest = tauFallFastDual * scaleFactors;
    tauFallSlowTest = tauFallSlowDual * scaleFactors;
    gVecs = compute_gabab_conductance(tVec, ipscStart, ampDual, tauRiseTest, ...
                                tauFallFastTest, tauFallSlowTest, weightDual);
    plot_conductance(tVec, gVecs, colorMap);
    tauEmpirical = siMs .* compute_time_constant(gVecs, 'DecayMethod', 'empirical');
    legend(create_labels_from_numbers(tauEmpirical, 'Prefix', 'tau = '));
    save_all_figtypes(fig, figPathBase, figTypes);
end

%% Plot as tauFallFastDual and tauFallSlowDual are varied with the same
%   relative value and weights
if plotVaryTauFall
    fig = set_figure_properties('AlwaysNew', true);
    figPathBase = fullfile(outFolder, 'gababipsc_dual_vary_taufall');
    scaleFactors = 0.1:0.1:1;
    tauFallFastTest = tauFallFastDual * scaleFactors;
    tauFallSlowTest = tauFallSlowDual * scaleFactors;
    gVecs = compute_gabab_conductance(tVec, ipscStart, ampDual, tauRiseDual, ...
                                tauFallFastTest, tauFallSlowTest, weightDual);
    plot_conductance(tVec, gVecs, colorMap);
    tauEmpirical = siMs .* compute_time_constant(gVecs, 'DecayMethod', 'empirical');
    legend(create_labels_from_numbers(tauEmpirical, 'Prefix', 'tau = '));
    save_all_figtypes(fig, figPathBase, figTypes);
end

%% Plot as tauFallFastDual and tauFallSlowDual are varied
%   while keeping the area under the curve constant
if plotVaryShapeOld
    fig = set_figure_properties('AlwaysNew', true);
    figPathBase = fullfile(outFolder, 'gababipsc_dual_vary_shape_old');
    tauFallFastTest = tauFallFastDual .* 2 .^ (-3:6);
    tauFallSlowInit = tauFallFastTest .* 5;
    tauFallSlowTest = compute_matching_tauFallSlow(aucDual, tauFallSlowInit, ...
                                    tVec, ipscStart, ampDual, tauRiseDual, ...
                                    tauFallFastTest, weightDual);
    gVecs = compute_gabab_conductance(tVec, ipscStart, ampDual, tauRiseDual, ...
                                tauFallFastTest, tauFallSlowTest, weightDual);
    tauEmpirical = siMs .* compute_time_constant(gVecs, 'DecayMethod', 'empirical');
    plot_conductance(tVec, gVecs, colorMap);
    legend(create_labels_from_numbers(tauEmpirical, 'Prefix', 'tau = '));
    save_all_figtypes(fig, figPathBase, figTypes);
end

%% Plot as all parameters are varied between Dual -> GAT3 shapes
%   while keeping the area under the curve constant
if plotVaryDualtoGAT3
    fig = set_figure_properties('AlwaysNew', true);
    figPathBase = fullfile(outFolder, 'gababipsc_vary_dual_to_gat3');
    scaleFactors = -0.8:0.2:1.4;
    ampTest = logscale(ampDual, ampGat3New, scaleFactors);
    tauRiseTest = logscale(tauRiseDual, tauRiseGat3, scaleFactors);
    tauFallFastTest = logscale(tauFallFastDual, tauFallFastGat3, scaleFactors);
    tauFallSlowInit = logscale(tauFallSlowDual, tauFallSlowGat3, scaleFactors);
    weightTest = logscale(weightDual, weightGat3, scaleFactors);
    tauFallSlowTest = compute_matching_tauFallSlow(aucDual, tauFallSlowInit, ...
                                    tVec, ipscStart, ampTest, tauRiseTest, ...
                                    tauFallFastTest, weightTest);
    gVecs = compute_gabab_conductance(tVec, ipscStart, ampTest, tauRiseTest, ...
                                tauFallFastTest, tauFallSlowTest, weightTest);
    tauEmpirical = siMs .* compute_time_constant(gVecs, 'DecayMethod', 'empirical');
    plot_conductance(tVec, gVecs, colorMap);
    legend(create_labels_from_numbers(tauEmpirical, 'Prefix', 'tau = '));
    save_all_figtypes(fig, figPathBase, figTypes);
end

%% Plot as all parameters are varied between Dual -> GAT3 -> GAT1 shapes
%   while keeping the area under the curve constant
if plotVaryDualtoGAT3toGAT1
    fig = set_figure_properties('AlwaysNew', true);
    figPathBase = fullfile(outFolder, 'gababipsc_vary_dual_to_gat3_to_gat1');
    scalesDualToGat3 = 0:0.2:1;
    ampTest = logscale(ampDual, ampGat3New, scalesDualToGat3);
    tauRiseTest = logscale(tauRiseDual, tauRiseGat3, scalesDualToGat3);
    tauFallFastTest = logscale(tauFallFastDual, tauFallFastGat3, scalesDualToGat3);
    tauFallSlowInit = logscale(tauFallSlowDual, tauFallSlowGat3, scalesDualToGat3);
    weightTest = logscale(weightDual, weightGat3, scalesDualToGat3);
    tauFallSlowTest = compute_matching_tauFallSlow(aucDual, tauFallSlowInit, ...
                                    tVec, ipscStart, ampTest, tauRiseTest, ...
                                    tauFallFastTest, weightTest);
    gVecsDualToGat3 = compute_gabab_conductance(tVec, ipscStart, ...
                                ampTest, tauRiseTest, ...
                                tauFallFastTest, tauFallSlowTest, weightTest);
    scalesGat3ToGat1 = 0.2:0.2:1;
    ampTest = logscale(ampGat3New, ampGat1New, scalesGat3ToGat1);
    tauRiseTest = logscale(tauRiseGat3, tauRiseGat1, scalesGat3ToGat1);
    tauFallFastTest = logscale(tauFallFastGat3, tauFallFastGat1, scalesGat3ToGat1);
    tauFallSlowInit = logscale(tauFallSlowGat3, tauFallSlowGat1, scalesGat3ToGat1);
    weightTest = logscale(weightGat3, weightGat1, scalesGat3ToGat1);
    tauFallSlowTest = compute_matching_tauFallSlow(aucDual, tauFallSlowInit, ...
                                    tVec, ipscStart, ampTest, tauRiseTest, ...
                                    tauFallFastTest, weightTest);
    gVecsGat3ToGat1 = compute_gabab_conductance(tVec, ipscStart, ...
                                ampTest, tauRiseTest, ...
                                tauFallFastTest, tauFallSlowTest, weightTest);
    gVecs = [gVecsDualToGat3, gVecsGat3ToGat1];
    tauEmpirical = siMs .* compute_time_constant(gVecs, 'DecayMethod', 'empirical');
    plot_conductance(tVec, gVecs, colorMap);
    legend(create_labels_from_numbers(tauEmpirical, 'Prefix', 'tau = '));
    save_all_figtypes(fig, figPathBase, figTypes);
end

%% Plot old figures
if plotOldFigures
    figs(2) = figure(56);
    pos(2,:) = defaultPosition;
    set(figs(2), 'Position', pos(2, :));
    tauFallFast = tauFallFastOrig;
    tauFallFast(3) = 1022.00;
    tauFallFast(4) = 2600.00;
    tauFallSlow = tauFallSlowOrig;
    tauFallSlow(3) = 273.40;
    tauFallSlow(4) = 65.80;
    weight = weightOrig;
    weight(3) = 0.225;
    weight(4) = 0.371;
    t0 = 3000;      % From "s.start = 3000" in build() of singleneuron4compgabab.hoc
    weight = 1;     % From "net = new NetCon(s, gababsyn, 0, 0, 1)" in build() of singleneuron4compgabab.hoc
    ipscStart = 0;      % Default of NMODL
    Rlast0 = 0;     % Default of NMODL
    Ninputs = 1;    % From "gababsyn.Ninputs = 1" in build() of singleneuron4compgabab.hoc
    Rdummy0 = 0;    % Default of NMODL
    RoffFast0 = 0;  % Default of NMODL
    RoffSlow0 = 0;  % Default of NMODL
    Ron0 = 0;       % Default of NMODL
    colorMap = {'k', 'b', 'r', 'g'};
    gVecs = computecond_weird(colorMap, tVec, ipscStart, ampOrig, ...
                                tauRiseOrig, tauFallFast, tauFallSlow, ...
                                weight, t0, actWeight, Rlast0, Ninputs, ...
                                Rdummy0, RoffFast0, RoffSlow0, Ron0);
    plot_conductance(tVec, gVecs, colorMap);
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
    title('GABAB IPSC conductances from optimizer4gabab and from dynamic clamp experiments 101210')
    xlim([0 4000])
    saveas(figs(2), fullfile(outFolder, 'gababipsc_optimizer_dclamp102210'), 'png')

    figs(3) = figure(57);
    pos(3,:) = defaultPosition;
    pos(3,1) = defaultPosition(1,1) + defaultPosition(1,3);
    set(figs(3), 'Position', pos(3,:));
    tauRise = [52.00; 52.00; 37.88; 39.9];            % (ms)
    tauFallFast = [90.10; 90.10; 962.4; 3010.0];      % (ms)
    tauFallSlow = [1073.20; 1073.20; 248.6; 65.8];    % (ms)
    weight = [0.952; 0.952; 0.247; 0.371];
    gVecs = compute_gabab_conductance(tVec, ipscStart, ampOrig, tauRise, ...
                                        tauFallFast, tauFallSlow, weight);
    plot_conductance(tVec, gVecs, colorMap);
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
    title('GABAB IPSC conductances from dynamic clamp experiments 091710')
    xlim([0 4000])
    saveas(figs(3), fullfile(outFolder, 'gababipsc_dclamp091710'), 'png')

    figs(4) = figure(58);
    pos(4,:) = defaultPosition;
    pos(4,1) = defaultPosition(1,1) - defaultPosition(1,3);
    pos(4,2) = defaultPosition(1,2) - defaultPosition(1,4);
    set(figs(4), 'Position', pos(4,:));
    amp = [8.888; 18.556; 12.752; 10.424];          % (nS)
    tauRise = [38.63; 45.74; 42.13; 40.18];           % (ms)
    tauFallFast = [1022; 127.5; 173.3; 217.5];        % (ms)
    tauFallSlow = [273.4; 1099; 1079; 1057];          % (ms)
    weight = [0.225; 0.913; 0.867; 0.827];
    gVecs = compute_gabab_conductance(tVec, ipscStart, amp, tauRise, ...
                                        tauFallFast, tauFallSlow, weight);
    plot_conductance(tVec, gVecs, colorMap);
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
    title('GABAB IPSC conductances from dynamic clamp experiments 102210')
    xlim([0 4000])
    saveas(figs(4), fullfile(outFolder, 'gababipsc_dclamp102210'), 'png')

    % TODO:
    % alignaxes(figs)

    figs(5) = figure(59);
    pos(5,:) = pos(4,:);
    set(figs(2), 'Position', pos(5,:));
    amp = [32.00; 48.00; 17.76; 12.64];                 % (nS)
    tauRise = [52.00; 52.00; 37.88; 39.9];              % (ms)
    tauFallFast = [90.10; 90.10; 962.4; 3010.0];        % (ms)
    tauFallSlow = [1073.20; 1073.20; 248.6; 65.8];      % (ms)
    weight = [0.952; 0.952; 0.247; 0.371];
    t0 = 0;         % start at time 0
    weight = 1;     %
    Ninputs = 1;    % 
    gVecs = computecond_m3ha(colorMap, tVec, ipscStart, amp, tauRise, ...
                    tauFallFast, tauFallSlow, weight, t0, actWeight, Ninputs);
    plot_conductance(tVec, gVecs, colorMap);
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
    title('GABAB IPSC conductances for gabab_m3ha.mod with 18 state variables')
    xlim([0 4000])
    saveas(figs(2), fullfile(outFolder, 'gababipsc_m3ha_network'), 'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gVecs = computecond_weird(colorMap, tVec, ipscStart, amp, tauRise, ...
                                    tauFallFast, tauFallSlow, weight, ...
                                    t0, actWeight, Rlast0, ...
                                    Ninputs, Rdummy0, RoffFast0, RoffSlow0, Ron0)
% Compute conductance curves that might have been implemented in optimizer4gabab
for k = 1:length(colorMap)
	temp1 = exp((ipscStart-t0)/(tauRise(k)));
	temp2 = exp((ipscStart-t0)/(tauFallFast(k)));
	temp3 = exp((ipscStart-t0)/(tauFallSlow(k)));
	Rlast1 = Rlast0 * ((1 - temp1)^8) * (weight(k)*temp2 + (1-weight(k))*temp3);
	
    Tlast1 = t0;
	act = (actWeight / Ninputs + Rlast0);
	Rlast2 = Rlast1 + act;
	Rdummy1 = Rdummy0 + Rlast1;
	RoffFast1 = RoffFast0 + act;
	RoffSlow1 = RoffSlow0 + act;
	Ron1 = Ron0 + act;

    phi = 1;
    Ron(k,:) = Ron1*exp(-phi*tVec/(tauRise(k)));
    RoffFast(k,:) = RoffFast1*exp(-phi*tVec/(tauFallFast(k)));
    RoffSlow(k,:) = RoffSlow1*exp(-phi*tVec/(tauFallSlow(k)));
    gVecs(:, k) = amp(k) * (1 - Ron(k,:)) .^ 8 ...
                .* (RoffFast(k,:)*weight(k) + RoffSlow(k,:)*(1-weight(k)));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gVecs = computecond_m3ha(colorMap, tVec, ipscStart, amp, tauRise, ...
                                    tauFallFast, tauFallSlow, weight, ...
                                    t0, actWeight, Ninputs)
% Compute conductance curves for new m3ha network using the Runge-Kutta method
for k = 1:length(colorMap)
    Tlast1 = t0;
	act = actWeight / Ninputs;
%{
    phi = 1;
    RFast0(k, :) = RFast0*exp(-phi*tVec/(tauRise(k)));
    RFast1(k, :) = RFast1*exp(-phi*tVec/(tauFallFast(k)));
    RFast2(k, :) = RSlow1*exp(-phi*tVec/(tauFallSlow(k)));

    gVecs(k, :) = amp(k) * (1 - Ron(k,:)) .^ 8 ...
                .* (RoffFast(k,:)*weight(k) + RoffSlow(k,:)*(1-weight(k)));  
%}
    [~, thisG] = ode45(odefun, tVec, ones(18, 1));
    gVecs(:, k) = thisG;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_conductance(tVec, gVecs, colorMap)
%% Plots conductance curves

plot_traces(tVec, gVecs, 'ColorMap', colorMap, ...
            'XLabel', 'Time (ms)', 'YLabel', 'Conductance (nS)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auc = compute_ipsc_area (tVec, ipscStart, amp, tauRise, ...
                                    tauFallFast, tauFallSlow, weight)

auc = trapz(compute_gabab_conductance(tVec, ipscStart, amp, tauRise, ...
                                    tauFallFast, tauFallSlow, weight));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tauFallSlow = compute_matching_tauFallSlow (gipscArea, ...
                            tauFallSlowInit, tVec, ipscStart, ...
                            amp, tauRise, tauFallFast, weight)

% Count the maximum number of conditions
% TODO for SHINSHIN: compute_maximum_numel.m
nConds = max([numel(tauFallSlowInit), numel(ipscStart), numel(amp), ...
                numel(tauRise), numel(tauFallFast), numel(weight)]);

% Match dimensions
[tauFallSlowInit, ipscStart, amp, tauRise, tauFallFast, weight] = ...
    argfun(@(x) match_dimensions(x, [nConds, 1]), ...
            tauFallSlowInit, ipscStart, amp, tauRise, tauFallFast, weight);

% Look for the matching tauFallSlow
tauFallSlow = ...
    array_fun(@(a, b, c, d, e, f) compute_one_matching_tauFallSlow(...
                                    gipscArea, tVec, a, b, c, d, e, f), ...
                tauFallSlowInit, ipscStart, amp, tauRise, tauFallFast, weight);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tauFallSlow = compute_one_matching_tauFallSlow (gipscArea, ...
                                        tVec, tauFallSlowInit, ipscStart, ...
                                        amp, tauRise, TFallFast, weight)

area_gipsc = ...
    @(x) compute_ipsc_area(tVec, ipscStart, amp, tauRise, TFallFast, x, weight);

sum_of_squared_error = @(x, y) (x - y) .* conj(x - y);

errorfun = @(x) sum_of_squared_error(area_gipsc(x), gipscArea);

tauFallSlow = fminsearch(errorfun, tauFallSlowInit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%