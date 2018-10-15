function [fitResults, fitObject, goodnessOfFit, algorithmInfo] = ...
                fit_pulse_response (tVec, yVec, pulseWidth, varargin)
%% Estimate short and long pulse response parameters from a double exponential fit to the rising/falling phase of a pulse response
% Usage: [fitResults, fitObject, goodnessOfFit, algorithmInfo] = ...
%               fit_pulse_response (tVec, yVec, pulseWidth, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       fitResults      - fitted parameters, including the fields: 
%                           eqnShortPulseResponse
%                           eqnLongPulseResponse
%                           ampSlowSpr  - amplitude of the slow component 
%                                               of the short pulse response
%                           tauSlow     - time constant of the slow component
%                           ampFastSpr  - amplitude of the fast component 
%                                               of the short pulse response
%                           tauFast     - time constant of the fast component
%                           ampSlowLpr  - amplitude of the slow component 
%                                               of the long pulse response
%                           ampFastLpr  - amplitude of the fast component 
%                                               of the long pulse response
%                           coeffSpr
%                           coeffLpr
%                           coeffSprStr
%                           coeffLprStr
%                           phaseName
%                           pulseWidth
%                           (see parse_fitobject.m for more)
%                       specified as a structure 
%       fitObject       - the cfit object returned by the fit function
%                       specified as a cfit object
%       goodnessOfFit   - the gof structure returned by the fit function
%                       specified as a gof structure 
%       algorithmInfo   - the output structure returned by fit_setup_2exp.m
%                       specified as a scalar structure
% Arguments:    
%       tVec        - time vector
%                   must be a numeric vector
%       yVec        - response vector
%                   must be a numeric vector
%       pulseWidth  - pulse width in the same units as time vector
%                   must be a nonnegative scalar
%       varargin    - 'PhaseName': the phase of the response
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto' - decide based on comparing absolute y values
%                       'rising'  - rising curve
%                       'falling' - falling curve
%                       'combined' - combined rising and falling curve
%                   default == 'auto'
%                   - 'AmplitudeEstimate': an estimated amplitude
%                   must be a numeric scalar
%                   default == []
%                   - 'MaxScalingFactor': maximum scaling factor for amplitude
%                   must be a nonnegative scalar
%                   default == 10
%                   - 'Tau1Init': estimated time constant for the 1st component
%                   must be a nonnegative scalar
%                   default == []
%                   - 'Tau2Init': estimated time constant for the 2nd component
%                   must be a nonnegative scalar
%                   default == []
%                   - 'Tau1Range': time constant range for the 1st component
%                   must be a numeric 2-element vector
%                   default == [0, Inf]
%                   - 'Tau2Range': time constant range for the 2nd component
%                   must be a numeric 2-element vector
%                   default == [0, Inf]
%
% Requires:
%       cd/fit_2exp.m
%
% Used by:    
%       cd/fit_and_estimate_passive_params.m

% File History:
% 2018-10-11 Moved code from find_passive_params.m
% 2018-10-11 Renamed coefficients to be more descriptive
% 2018-10-14 Added the combined phase
% 

%% Hard-coded parameters
validPhaseNames = {'rising', 'falling', 'combined', 'auto'};
nSigFig = 2;    % maximum number of significant figures in the equation string
eqFormCombinedTemplate = ['( a*(1-exp(-x/b)) + ', ...
                            'c*(1-exp(-x/d)) ) * (x <= pw) + ', ...
                            '( a*(1-exp(-pw/b))*exp(-(x - pw)/b) + ', ...
                            'c*(1-exp(-pw/d))*exp(-(x - pw)/d) ) * (x > pw)'];

%% Default values for optional arguments
phaseNameDefault = 'auto';          % set later
amplitudeEstimateDefault = [];
maxScalingFactorDefault = 10;
tau1InitDefault = [];
tau2InitDefault = [];
tau1RangeDefault = [0, Inf];
tau2RangeDefault = [0, Inf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'tVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'yVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'pulseWidth', ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PhaseName', phaseNameDefault, ...
    @(x) any(validatestring(x, validPhaseNames)));
addParameter(iP, 'AmplitudeEstimate', amplitudeEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'MaxScalingFactor', maxScalingFactorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Tau1Init', tau1InitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Tau2Init', tau2InitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Tau1Range', tau1RangeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'Tau2Range', tau2RangeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));

% Read from the Input Parser
parse(iP, tVec, yVec, pulseWidth, varargin{:});
phaseName = validatestring(iP.Results.PhaseName, validPhaseNames);
amplitudeEstimate = iP.Results.AmplitudeEstimate;
maxScalingFactor = iP.Results.MaxScalingFactor;
tau1Init = iP.Results.Tau1Init;
tau2Init = iP.Results.Tau2Init;
tau1Range = iP.Results.Tau1Range;
tau2Range = iP.Results.Tau2Range;

% Check relationships between arguments
if length(tVec) ~= length(yVec)
    error('tVec and yVec must have the same length!!');
end

%% Preparation
% Estimate a direction if not provided
if strcmpi(phaseName, 'auto')
    % Extract the absolute value of the first and last y value
    yFirstAbs = abs(yVec(1));
    yLastAbs = abs(yVec(end));

    % Extract the y value that has the largest absolute value
    yMaxAbs = max(abs(yVec));

    % Check if its a combined phase by checking whether
    %   the average of the first and last y values are less than
    %   1/4 of the maximum absolute y value
    % Otherwise, check whether the first y value is closer or farther away 
    %   from zero than the last y value
    if mean([yFirstAbs, yLastAbs]) < (1/4) * yMaxAbs
        phaseName = 'combined';
    elseif yFirstAbs <= yLastAbs
        phaseName = 'rising';
    else
        phaseName = 'falling';
    end
end

% Decide on the direction and any custom equation form
switch phaseName
    case 'combined'
        direction = 'custom';
        eqFormCustom = ...
            replace(eqFormCombinedTemplate, 'pw', num2str(pulseWidth, 2));
    case {'rising', 'falling'}
        direction = phaseName;
        eqFormCustom = '';
    otherwise
        error('phaseName unrecognized!!');
end


%% Fit to a double exponential curve
[fitResults, fitObject, goodnessOfFit, algorithmInfo] = ...
    fit_2exp(tVec, yVec, ...
            'Direction', direction, ...
            'AmplitudeEstimate', amplitudeEstimate, ...
            'MaxScalingFactor', maxScalingFactor, ...
            'Tau1Init', tau1Init, ...
            'Tau2Init', tau2Init, ...
            'Tau1Range', tau1Range, ...
            'Tau2Range', tau2Range, ...
            'EquationForm', eqFormCustom);

% Extract from fitResults
fitFormula = fitResults.fitFormula;
coeffNames = fitResults.coeffNames;
coeffValues = fitResults.coeffValues;

%% Estimate coefficients
% Rename the short pulse response coefficients 
%   to be consistent with Johnston & Wu, p. 94
%   Note: 1. Rather than 0 & 1, use slow & fast for clarity
%         2. Force the first component to have the larger time constant
if coeffValues(2) >= coeffValues(4)
    ampSlowSpr = coeffValues(1);
    tauSlow = coeffValues(2);
    ampFastSpr = coeffValues(3);
    tauFast = coeffValues(4);
else
    ampSlowSpr = coeffValues(3);
    tauSlow = coeffValues(4);
    ampFastSpr = coeffValues(1);
    tauFast = coeffValues(2);
end

% Put the short pulse response coefficients together
coeffSpr = [ampSlowSpr; tauSlow; ampFastSpr; tauFast];

% Compute corresponding long pulse response coefficients
[ampSlowLpr, ampFastLpr] = ...
    coeffSpr2coeffLpr(ampSlowSpr, tauSlow, ampFastSpr, tauFast, ...
                        phaseName, pulseWidth);

% Put the long pulse response coefficients together
coeffLpr = [ampSlowLpr; tauSlow; ampFastLpr; tauFast];

%% Generate equations
% Generate strings for fitted coefficients
coeffSprStr = arrayfun(@(x) num2str(x, nSigFig), coeffSpr, ...
                                'UniformOutput', false);
coeffLprStr = arrayfun(@(x) num2str(x, nSigFig), coeffLpr, ...
                                'UniformOutput', false);

% Generate strings for the fitted equations
eqnShortPulseResponse = replace(fitFormula, coeffNames, coeffSprStr);
eqnLongPulseResponse = replace(fitFormula, coeffNames, coeffLprStr);

%% Place outputs in fitResults structure
fitResults.eqnShortPulseResponse = eqnShortPulseResponse;
fitResults.eqnLongPulseResponse = eqnLongPulseResponse;
fitResults.ampSlowSpr = ampSlowSpr;
fitResults.tauSlow = tauSlow;
fitResults.ampFastSpr = ampFastSpr;
fitResults.tauFast = tauFast;
fitResults.ampSlowLpr = ampSlowLpr;
fitResults.ampFastLpr = ampFastLpr;
fitResults.coeffSpr = coeffSpr;
fitResults.coeffLpr = coeffLpr;
fitResults.coeffSprStr = coeffSprStr;
fitResults.coeffLprStr = coeffLprStr;
fitResults.phaseName = phaseName;
fitResults.pulseWidth = pulseWidth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ampSlowLpr, ampFastLpr] = ...
                coeffSpr2coeffLpr (ampSlowSpr, tauSlow, ampFastSpr, ...
                                    tauFast, phaseName, pulseWidth)
%% Estimate the long pulse response coefficients from the 
%   short pulse response coefficients
% See equation 4.5.56 on p. 95 of Johnston % Wu

% Depends on phase
switch phaseName
    case 'falling'
        % The short pulse coefficients are "partially charged"
        %   from the long pulse coefficients
        %   Thus, one can back-calculate the long pulse coefficients
        ampSlowLpr = ampSlowSpr / (1 - exp(-pulseWidth/tauSlow));
        ampFastLpr = ampFastSpr / (1 - exp(-pulseWidth/tauFast));        
    case {'rising', 'combined'}
        % The long pulse coefficients should be the same as 
        %   the short pulse coefficients
        ampSlowLpr = ampSlowSpr;
        ampFastLpr = ampFastSpr;
    otherwise
        error('phaseName unrecognized!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if fitObject.b >= fitObject.d
    ampSlowSpr = fitObject.a;
    tauSlow = fitObject.b;
    ampFastSpr = fitObject.c;
    tauFast = fitObject.d;
else
    ampSlowSpr = fitObject.c;
    tauSlow = fitObject.d;
    ampFastSpr = fitObject.a;
    tauFast = fitObject.b;
end

if strcmp(phaseName, 'falling')
    [ampSlowLpr, ampFastLpr] = coeffSpr2coeffLpr(ampSlowSpr, tauSlow, ampFastSpr, tauFast, pulseWidth);
elseif strcmp(phaseName, 'rising')
    ampSlowLpr = ampSlowSpr;
    ampFastLpr = ampFastSpr;
else
    ampSlowLpr = NaN;
    ampFastLpr = NaN;
end

eqnShortPulseResponse = strrep(equationForm, 'a', num2str(ampSlowSpr, nSigFig));
eqnShortPulseResponse = strrep(eqnShortPulseResponse, 'b', num2str(tauSlow, nSigFig));
eqnShortPulseResponse = strrep(eqnShortPulseResponse, 'c', num2str(ampFastSpr, nSigFig));
eqnShortPulseResponse = strrep(eqnShortPulseResponse, 'd', num2str(tauFast, nSigFig));
eqnLongPulseResponse = strrep(equationForm, 'a', num2str(ampSlowLpr, nSigFig));
eqnLongPulseResponse = strrep(eqnLongPulseResponse, 'b', num2str(tauSlow, nSigFig));
eqnLongPulseResponse = strrep(eqnLongPulseResponse, 'c', num2str(ampFastLpr, nSigFig));
eqnLongPulseResponse = strrep(eqnLongPulseResponse, 'd', num2str(tauFast, nSigFig));

nCoeffs = fitResults.nCoeffs;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%