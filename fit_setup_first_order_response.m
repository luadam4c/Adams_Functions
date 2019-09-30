function [eqForm, aFittype, coeffInit, coeffLower, coeffUpper] = ...
                fit_setup_first_order_response (varargin)
%% Constructs a fittype object and set initial conditions and bounds for a first order response equation form
% Usage: [eqForm, aFittype, coeffInit, coeffLower, coeffUpper] = ...
%               fit_setup_first_order_response (varargin)
% Explanation:
%   This function is for fitting equations of the form
%       a*(1-exp(-(x - d)/b))*(x > d & x <= d + c) + 
%       a*(1-exp(-c/b))*exp(-(x - d - c)/b)*(x > d + c)
%
% Example(s):
%       fit_setup_first_order_response;
%       fit_setup_first_order_response('MaxScalingFactor', 10);
%       fit_setup_first_order_response('TotalDuration', 100);
%
% Outputs:
%       eqForm      - equation form used
%                   specified as a character vector
%       aFittype    - fittype object constructed
%                   specified as a fittype object
%       coeffInit   - initial values for the coefficients
%                   specified as a numeric column vector
%       coeffLower  - lower bounds for the coefficients
%                   specified as a numeric column vector
%       coeffUpper  - upper bounds for the coefficients
%                   specified as a numeric column vector
%
% Arguments:    
%       varargin    - 'TotalDuration': total duration for the response
%                   must be a numeric scalar
%                   default == Inf
%                   - 'AmplitudeEstimate': an estimated amplitude
%                   must be a numeric scalar
%                   default == 1
%                   - 'MaxScalingFactor': maximum scaling factor for amplitude
%                   must be a nonnegative scalar
%                   default == 10
%                   - 'TauEstimate': estimated time constant
%                   must be a nonnegative scalar
%                   default == TotalDuration / 3
%                   - 'TauRange': time constant range
%                   must be a numeric 2-element vector
%                   default == [0, Inf]
%                   - 'DurationEstimate': estimated stimulus duration
%                   must be a nonnegative scalar
%                   default == TotalDuration / 3
%                   - 'DelayEstimate': estimated delay
%                   must be a numeric 2-element vector
%                   default == TotalDuration / 3
%                   - 'EquationForm': equation form if the direction is 'custom'
%                   must be a character vector
%                   default == ''
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/fit_first_order_response.m TODO

% File History:
% 2019-09-09 Modified from fit_setup_2exp.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
totalDurationDefault = Inf;
amplitudeEstimateDefault = 1;
maxScalingFactorDefault = 10;
tauEstimateDefault = [];        % set later
tauRangeDefault = [0, Inf];
durationEstimateDefault = [];   % set later
delayEstimateDefault = [];      % set later
eqFormDefault = ['a*(1-exp(-(x - d)/b))*(x > d & x <= d + c) + ', ...
                'a*(1-exp(-c/b))*exp(-(x - d - c)/b)*(x > d + c)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TotalDuration', totalDurationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'AmplitudeEstimate', amplitudeEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MaxScalingFactor', maxScalingFactorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'TauEstimate', tauEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'TauRange', tauRangeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'DurationEstimate', durationEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'DelayEstimate', delayEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'EquationForm', eqFormDefault, ...
    @(x) validateattributes(x, {'char'}, {'2d'}));

% Read from the Input Parser
parse(iP, varargin{:});
totalDuration = iP.Results.TotalDuration;
amplitudeEstimate = iP.Results.AmplitudeEstimate;
maxScalingFactor = iP.Results.MaxScalingFactor;
tauEstimate = iP.Results.TauEstimate;
tauRange = iP.Results.TauRange;
durationEstimate = iP.Results.DurationEstimate;
delayEstimate = iP.Results.DelayEstimate;
eqForm = iP.Results.EquationForm;

%% Preparation
% Compute the default time constant estimate
if isempty(tauEstimate)
    if ~isinf(totalDuration)
        tauEstimate = totalDuration / 3;
    else
        tauEstimate = 1;
    end
end

% Compute the default duration estimate
if isempty(durationEstimate)
    if ~isinf(totalDuration)
        durationEstimate = totalDuration / 3;
    else
        durationEstimate = 1;
    end
end

% Compute the default delay estimate
if isempty(delayEstimate)
    if ~isinf(totalDuration)
        delayEstimate = totalDuration / 3;
    else
        delayEstimate = 1;
    end
end

%% Do the job
% Construct a fittype object based on an equation form
aFittype = fittype(eqForm);

% Get all the coefficient names
%   Note: This is most likely: {'a'; 'b'; 'c'; 'd'}, 
%           but not necessarily so in this order
coeffNames = coeffnames(aFittype);

% Count the number of coefficients
nCoeffs = numel(coeffNames);

% Set up initial conditions and boundary conditions for the coefficients
coeffInit = zeros(nCoeffs, 1);      % initial guesses for the coefficients
coeffLower = zeros(nCoeffs, 1);     % lower bounds for the coefficients
coeffUpper = zeros(nCoeffs, 1);     % upper bounds for the coefficients
for iCoeff = 1:nCoeffs
    % Get the current coefficient name
    coeffName = coeffNames{iCoeff};

    % Decide on the coefficient initial values and bounds
    switch coeffName
    case 'a'                        % amplitude
        coeffInit(iCoeff) = amplitudeEstimate;
        if amplitudeEstimate > 0
            coeffLower(iCoeff) = 0;
            coeffUpper(iCoeff) = maxScalingFactor * amplitudeEstimate;
        else
            coeffLower(iCoeff) = maxScalingFactor * amplitudeEstimate;
            coeffUpper(iCoeff) = 0;
        end
    case 'b'                        % time constant
        coeffInit(iCoeff) = tauEstimate;
        coeffLower(iCoeff) = tauRange(1);
        coeffUpper(iCoeff) = tauRange(2);
    case 'c'                        % pulse width
        coeffInit(iCoeff) = durationEstimate;
        coeffLower(iCoeff) = 0;
        coeffUpper(iCoeff) = totalDuration;
    case 'd'                        % delay
        coeffInit(iCoeff) = delayEstimate;
        coeffLower(iCoeff) = 0;
        coeffUpper(iCoeff) = totalDuration;
    otherwise
        error('coeffName is unrecognized!!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%