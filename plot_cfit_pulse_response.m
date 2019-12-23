function h = plot_cfit_pulse_response (xVec, yVec, varargin)
%% Plots data along with the fitted curve
% Usage: h = plot_cfit_pulse_response (xVec, yVec, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - figure handle for the created figure
%                   specified as a figure handle
% Arguments:
%       xVec        - x vector used for the fit
%                   must be a numeric vector
%       yVec        - y vector used for the fit
%                   must be a numeric vector
%       varargin    - 'FitObject': the cfit object returned by the 
%                                   fit() function
%                   must be a cfit object
%                   default == returned by fit_and_estimate_passive_params.m
%                   - 'FitResults': a structure containing the fitting results
%                   must be a structure containing at least the fields:
%                       eqnShortPulseResponse
%                       eqnLongPulseResponse
%                       tauSlow
%                       tauFast
%                       ampSlowLpr
%                       ampFastLpr
%                       phaseName
%                       pulseWidth
%                   default == returned by fit_and_estimate_passive_params.m
%                   - 'GoodnessOfFit': a structure containing 
%                                       goodness of fit measures
%                   must be a structure containing at least the fields:
%                       rmse
%                   default == returned by fit_and_estimate_passive_params.m
%                   - 'PassiveParams': a structure containing the 
%                                       estimated passive parameters
%                   must be a structure containing at least the fields:
%                       responseAmplitude
%                       Rinput
%                   default == returned by fit_and_estimate_passive_params.m
%                   - 'PulseWidth': pulse width in the same units as time vector
%                   must be a nonnegative scalar
%                   default == must be provided if FitObject is not provided
%                   - 'PulseAmplitude': pulse amplitude in the same units 
%                                       as x vector
%                   must be a numeric scalar
%                   default == must be provided if FitObject is not provided
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [min(timeVec), max(timeVec)]
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'northeast'
%                   - 'Component1Color': color for component 1
%                   must be a character vector or 3-element numeric array
%                   default == rgb('Turquoise')
%                   - 'Component2Color': color for component 2
%                   must be a character vector or 3-element numeric array
%                   default == rgb('DarkGreen')
%
% Requires:
%       cd/argfun.m
%       cd/compute_axis_limits.m
%       cd/create_latex_string.m
%       cd/fit_and_estimate_passive_params.m
%       cd/force_column_vector.m
%       cd/islegendlocation.m
%
% Used by:    
%       cd/find_passive_params.m
%       cd/m3ha_plot_figure03.m

% File History:
% 2018-10-12 Modified from find_passive_params.m
% 2018-10-12 Now uses phaseName stored in fit_pulse_response.m
% 2018-10-15 Added display of root-mean-square error
% 

%% Hard-coded parameters
yCoverage = 90;                 % coverage of y axis (%)
nSigFig = 3;                    % number of significant figures for Rinput
sprColor = [0.8594, 0.0781, 0.2344];    % rgb('Crimson');
lprColor = [0.7266, 0.3320, 0.8242];    % rgb('MediumOrchid');
rmseColor = [0, 0, 0.5000];             % rgb('Navy');
asymptoteColor = [0, 0, 0.5000];        % rgb('Navy');
RinColor = [0.2930, 0, 0.5078];         % rgb('Indigo');
lineWidth = 2;

%% Default values for optional arguments
fitObjectDefault = [];          % set later
fitResultsDefault = [];         % set later
goodnessOfFitDefault = [];      % set later
passiveParamsDefault = [];      % set later
pulseWidthDefault = [];         % set later
pulseAmplitudeDefault = [];     % set later
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later
legendLocationDefault = 'northeast';
component1ColorDefault = [0.2500, 0.8750, 0.8125];  % rgb('Turquoise');
component2ColorDefault = [0, 0.3906, 0];            % rgb('DarkGreen');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'xVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'yVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FitObject', fitObjectDefault, ...
    @(x) validateattributes(x, {'cfit'}, {'2d'}));
addParameter(iP, 'FitResults', fitResultsDefault, ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));
addParameter(iP, 'GoodnessOfFit', goodnessOfFitDefault, ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));
addParameter(iP, 'PassiveParams', passiveParamsDefault, ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));
addParameter(iP, 'PulseWidth', pulseWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'PulseAmplitude', pulseAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'Component1Color', component1ColorDefault);
addParameter(iP, 'Component2Color', component2ColorDefault);

% Read from the Input Parser
parse(iP, xVec, yVec, varargin{:});
fitObject = iP.Results.FitObject;
fitResults = iP.Results.FitResults;
goodnessOfFit = iP.Results.GoodnessOfFit;
passiveParams = iP.Results.PassiveParams;
pulseWidth = iP.Results.PulseWidth;
pulseAmplitude = iP.Results.PulseAmplitude;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
component1Color = iP.Results.Component1Color;
component2Color = iP.Results.Component2Color;

% Check relationships between arguments
% TODO

%% Extract from arguments
% Fit data if not yet done
if isempty(fitObject)
    % Make sure pulse parameters are provided
    if isempty(pulseWidth)
        fprintf('PulseWidth must be provided if FitObject not provided!!\n');
        h = [];
        return
    elseif isempty(pulseAmplitude)
        fprintf(['PulseAmplitude must be provided if ', ...
                'FitObject not provided!!\n']);
        h = [];
        return
    end

    % Provide warnings if fitResults and passiveParams already provided
    if ~isempty(fitResults)
        fprintf(['Warning: Provided fitResults will ', ...
                    'be overwritten by new fitting!\n']);
    end
    if ~isempty(passiveParams)
        fprintf(['Warning: Provided passiveParams will ', ...
                    'be overwritten by new fitting!\n']);
    end
    if ~isempty(goodnessOfFit)
        fprintf(['Warning: Provided goodnessOfFit will ', ...
                    'be overwritten by new fitting!\n']);
    end

    % Fit the data
    [passiveParams, fitResults, fitObject, goodnessOfFit] = ...
        fit_and_estimate_passive_params(xVec, yVec, pulseWidth, pulseAmplitude);
else
    if isempty(fitResults)
        fprintf(['FitResults must be provided if ', ...
                'FitObject provided!!\n']);
        h = [];
        return
    elseif isempty(passiveParams)
        fprintf(['PassiveParams must be provided if ', ...
                'FitObject provided!!\n']);
        h = [];
        return
    elseif isempty(goodnessOfFit)
        fprintf(['GoodnessOfFit must be provided if ', ...
                'FitObject provided!!\n']);
        h = [];
        return
    end
end

% Extract from fitResults
eqnS = fitResults.eqnShortPulseResponse;
eqnL = fitResults.eqnLongPulseResponse;
tau1 = fitResults.tauSlow;
tau2 = fitResults.tauFast;
ampRising1 = fitResults.ampSlowLpr;
ampRising2 = fitResults.ampFastLpr;
phaseName = fitResults.phaseName;
pulseWidth = fitResults.pulseWidth;

% Extract from passiveParams
responseAmplitude = passiveParams.responseAmplitude;
Rinput = passiveParams.Rinput;

% Extract from goodnessOfFit
rmse = goodnessOfFit.rmse;

%% Preparation
% Hold on for multiple plots
hold on

% Find the minimum and maximum time of interest
minT = min(xVec);
maxT = max(xVec);

% Get the sampling interval
siMs = xVec(2) - xVec(1);

% Compute the number of samples for the rising phase
nSamplesRising = round(pulseWidth / siMs);

% Find the minimum and maximum y values of interest
minY = min([responseAmplitude, 0]);
maxY = max([responseAmplitude, 0]);

% Set default x limits as the entire time range
if isempty(xLimits)
    xLimits = [minT, maxT];
end

% Set default y limits using the responseAmplitude
if isempty(yLimits)
    % Compute the y limits
    yLimits = compute_axis_limits([minY, maxY], 'y', 'Coverage', yCoverage);
end

% Break up the equation for the combined phase
% TODO: Make a function break_equation.m
if strcmp(phaseName, 'combined')
    % Find the index of the first ' + '
    idxFirstConnection = strfind(eqnS, ' + ');

    % Get the first part of the equation
    eqnSPart1 = eqnS(1:(idxFirstConnection-1));

    % Get the rest of the equation
    eqnSRest = eqnS((idxFirstConnection+1):end);

    % Find the index of the second ' + '
    idxSecondConnection = strfind(eqnSRest, ' + ');

    % Get the second part of the equation
    eqnSPart2 = eqnSRest(1:(idxSecondConnection-1));

    % Get the rest of the equation
    eqnSPart3 = eqnSRest((idxSecondConnection+1):end);

else
    eqnSPart1 = '';
    eqnSPart2 = '';
    eqnSPart3 = '';
end

% Generate latex strings for each equation
[eqnSLatex, eqnLLatex, eqnSPart1Latex, eqnSPart2Latex, eqnSPart3Latex] = ...
    argfun(@create_latex_string, eqnS, eqnL, eqnSPart1, eqnSPart2, eqnSPart3);

%% Plotting
% Plot the fitted curve with the data
hObjects = plot(fitObject, xVec, yVec);

% The fitted curve is the second object
fittedCurve = hObjects(2);

% Change the linewidth of the fitted curve
set(fittedCurve, 'LineWidth', lineWidth);

% Remove the default legend
legend('off');

% Create a new legend
if ~strcmpi(legendLocation, 'suppress')
    legend('data', 'fitted', 'Location', legendLocation);
end

% Plot each component separately
switch phaseName
case 'rising'
    % Generate a column time vector for the rising phase
    xVecRising = transpose(linspace(0, maxT, nSamplesRising));

    % Generate response vectors for the rising phase
    yVecRisingComp1 = ampRising1 * (1 - exp(-xVecRising/tau1));
    yVecRisingComp2 = ampRising2 * (1 - exp(-xVecRising/tau2));

    % Plot the components
    plot(xVecRising, yVecRisingComp1, 'Color', component1Color, ...
            'LineStyle', '--', 'DisplayName', 'Comp1', ...
            'LineWidth', lineWidth);
    plot(xVecRising, yVecRisingComp2, 'Color', component2Color, ...
            'LineStyle', '--', 'DisplayName', 'Comp2', ...
            'LineWidth', lineWidth);
case {'falling', 'combined'}
    % Generate time vectors for the falling and combined phases
    if strcmpi(phaseName, 'falling')
        % Generate a column time vector for the rising phase
        xVecRising = transpose(linspace(-pulseWidth, 0, nSamplesRising));

        % Compute the number of samples for the falling phase
        nSamplesFalling = round((maxT - minT) / siMs);

        % Generate a column time vector for the falling phase
        xVecFalling = transpose(linspace(0, maxT, nSamplesFalling));

        % Combine the rising and falling phases
        xVecCombined = [xVecRising; xVecFalling];
    else
        % Compute the number of samples for the combined phases
        nSamplesCombined = round((maxT - minT) / siMs);

        % Generate a column time vector for the combined phases
        xVecCombined = transpose(linspace(0, maxT, nSamplesCombined));

        % Extract the x vector for the rising phase
        xVecRising = xVecCombined(1:nSamplesRising);

        % Extract the x vector for the falling phase
        xVecFalling = xVecCombined(nSamplesRising+1:end);
    end

    % Generate response vectors for the rising phase
    xVecRisingShifted = xVecRising - xVecRising(1);
    yVecRisingComp1 = ampRising1 * (1 - exp(-xVecRisingShifted/tau1));
    yVecRisingComp2 = ampRising2 * (1 - exp(-xVecRisingShifted/tau2));

    % Record the final values for each component
    ampFalling1 = ampRising1 * (1 - exp(-pulseWidth/tau1));
    ampFalling2 = ampRising2 * (1 - exp(-pulseWidth/tau2));

    % Generate response vectors for the falling phase
    xVecFallingShifted = xVecFalling - xVecFalling(1);
    yVecFallingComp1 = ampFalling1 * exp(-xVecFallingShifted/tau1);
    yVecFallingComp2 = ampFalling2 * exp(-xVecFallingShifted/tau2);

    % Combine the rising and falling phases
    yVecComp1 = [yVecRisingComp1; yVecFallingComp1];
    yVecComp2 = [yVecRisingComp2; yVecFallingComp2];

    % Plot the components
    plot(xVecCombined, yVecComp1, 'Color', component1Color, ...
            'LineStyle', '--', 'DisplayName', 'Comp1', ...
            'LineWidth', lineWidth);
    plot(xVecCombined, yVecComp2, 'Color', component2Color, ...
            'LineStyle', '--', 'DisplayName', 'Comp2', ...
            'LineWidth', lineWidth);

    % Set new x-axis limits
    xLimits = [min(xLimits(1), min(xVecRising)), ...
                max(xLimits(2), max(xVecFalling))];
otherwise
    error('phaseName unrecognized!!');
end

% Set time axis limits
if ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end

% Set y axis limits
if ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

% Plot a line for the asymptotic response
line(xLimits, responseAmplitude * [1, 1], ...
        'Color', asymptoteColor, 'LineStyle', '--', 'LineWidth', 0.5, ...
        'DisplayName', 'Amp');

% Define starting x and y positions for texts
xpos = 1/30;
switch phaseName
case 'rising'
    % 3 lines are needed
    ypos = 5/30;
case 'falling'
    % 4 lines are needed
    ypos = 7/30;
case 'combined'
    % 5 lines are needed
    ypos = 9/30;
otherwise
end

% Show the root-mean-square error
text(xpos, ypos, ['root-mean-square error = ', num2str(rmse), ' mV'], ...
    'FontSize', 8, 'Color', rmseColor, ...
    'Units', 'normalized');

% Show the short-pulse response equation for all phases
%   but break up the equation into two for the combined phase
if strcmp(phaseName, 'rising') || strcmp(phaseName, 'falling')
    ypos = ypos - (1/15);
    text(xpos, ypos, eqnSLatex, ...
        'FontSize', 8, 'Color', sprColor, ...
        'Interpreter', 'latex', 'Units', 'normalized');
elseif strcmp(phaseName, 'combined')
    % Show the first part of the equation
    ypos = ypos - (1/15);
    text(xpos, ypos, eqnSPart1Latex, ...
        'FontSize', 8, 'Color', sprColor, ...
        'Interpreter', 'latex', 'Units', 'normalized');

    % Show the second part of the equation
    ypos = ypos - (1/15);
    text(xpos, ypos, eqnSPart2Latex, ...
        'FontSize', 8, 'Color', sprColor, ...
        'Interpreter', 'latex', 'Units', 'normalized');

    % Show the third part of the equation
    ypos = ypos - (1/15);
    text(xpos, ypos, eqnSPart3Latex, ...
        'FontSize', 8, 'Color', sprColor, ...
        'Interpreter', 'latex', 'Units', 'normalized');
end

% Show the long-pulse response equation for the falling phase only
if strcmp(phaseName, 'falling')
    ypos = ypos - (1/15);
    text(xpos, ypos, eqnLLatex, ...
        'FontSize', 8, 'Color', lprColor, ...
        'Interpreter', 'latex', 'Units', 'normalized');
end

% Show Rinput
ypos = ypos - (1/15);
text(xpos, ypos, ...
    ['Rin = ', num2str(Rinput, nSigFig), ' MOhm'], ...
    'FontSize', 8, 'Color', RinColor, ...
    'Units', 'normalized');

% Add axes labels
xlabel('Shifted time (ms)')
ylabel('Shifted voltage (mV)')

% TODO
h = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%function plot_cfit_pulse_response (phaseName, xLimits, fitObject, xvec, yvec, params, eqnS, eqnL, tau1, tau2, CL1, CL2, pulseWidth)

xVec2 = floor(-timeToPad):1:400;
yVec2 = -ampRising1*exp(-(xVec2+timeToPad)/tau1) - ampRising2*exp(-(xVec2+timeToPad)/tau2);

% Adjust the axes
ylim([(responseAmplitude - 0.5) 0.5]);        % Add 0.5 mV above and below
if strcmpi(phaseName, 'rising')
    xlim(xLimits - xLimits(1));
elseif strcmpi(phaseName, 'falling')
    xlim([-timeToPad xLimits(2)]);
end

% Get axes limits
ax = gca;
xLimits = get(ax, 'Xlim');
yLimits = get(ax, 'Ylim');

yVec2 = -ampRising1*exp(-xVec2Shifted/tau1) - ampRising2*exp(-xVec2Shifted/tau2);

% Get axes ranges
xRange = xLimits(2) - xLimits(1);
yRange = yLimits(2) - yLimits(1);

xpos = xLimits(1) + (1/30) * xRange;
ypos = yLimits(1) + (5/30) * yRange;
text(xpos, ypos, eqnS, 'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (1/15) * yRange, eqnL, ...
    'FontSize', 8, 'Color', phaseColor);
text(xpos, ypos - (2/15) * yRange, ...
    ['Rin = ', num2str(Rinput, nSigFig), ' MOhm'], ...
    'FontSize', 8, 'Color', phaseColor);

% Plot the long pulse response that corresponds to the short pulse response
%   for the falling phase
if strcmpi(phaseName, 'falling')
    % Get the sampling interval
    siMs = xVec(2) - xVec(1);

    % Get the maximum time
    maxTime = max(xVec);

    % Get the number of samples
    nSamples = length(xVec);

    % Find time before the short pulse response 
    %       the corresponding long pulse response would begin
    %   CSn = CLn * (1 - exp(-w/taun)) = CLn * exp(-t/taun))
    %   => t = -log(1 - exp(-w/taun)) * taun
    timeToPadComp1 = -log(1 - exp(-pulseWidth/tau1)) * tau1;
    timeToPadComp2 = -log(1 - exp(-pulseWidth/tau2)) * tau2;
    % fprintf('timeToPadComp1 == %g\n', timeToPadComp1);
    % fprintf('timeToPadComp2 == %g\n', timeToPadComp2);

    % Get the amount of padded time needed
    timeToPad = max(timeToPadComp1, timeToPadComp2);

    % Compute the number of samples to padd
    samplesToPad = round(timeToPad / siMs);

    % Compute the new number of samples
    nSamples2 = nSamples + samplesToPad;

    % Generate a new time vector
    xVec2 = linspace(-timeToPad, maxTime, nSamples2);

    % Generate a shifted time vector
    xVec2Shifted = xVec2 + timeToPad;

    % Generate a new response vector
    yVec2 = ampRising1*exp(-xVec2Shifted/tau1) + ampRising2*exp(-xVec2Shifted/tau2);

    % Plot the corresponding long pulse response
    plot(xVec2, yVec2, 'Color', phaseColor, 'DisplayName', 'LPR');

    % Set new x-axis limits
    xLimits = [min(xLimits(1), min(xVec2)), ...
                max(xLimits(2), max(xVec))];
end

legend('location', 'northeast');

% Decide on a phase-dependent color
if strcmp(phaseName, 'rising')
    phaseColor = risingColor;
elseif strcmp(phaseName, 'falling')
    phaseColor = fallingColor;
end

risingColor = 'r';
fallingColor = 'm';

if strcmpi(phaseName, 'rising')
    legend('data', 'fitted', 'Location', legendLocation);
elseif strcmpi(phaseName, 'falling')
    legend('data', 'fitted', 'Location', legendLocation);
end

% The x vector is for the rising phase
xVecRising = xVec;

% The x vector is for the falling phase
%   Note: make sure it is a column vector
xVecFalling = force_column_vector(xVec);

% The x vector is already for the combined phases
%   Note: make sure it is a column vector
xVecCombined = force_column_vector(xVec);

ypos = 5/30;

% Show the short-pulse response equation
switch phaseName
case {'rising', 'combined'}
    yposNow = ypos - (1/15);
case 'falling'
    yposNow = ypos;
otherwise
end

yLimits = compute_ylimits(minY, maxY, 'Coverage', yCoverage);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
