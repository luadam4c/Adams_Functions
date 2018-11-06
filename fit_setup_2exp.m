function [eqForm, aFittype, coeffInit, coeffLower, coeffUpper] = ...
                fit_setup_2exp (direction, varargin)
%% Constructs a fittype object and set initial conditions and bounds for a double exponential equation form
% Usage: [eqForm, aFittype, coeffInit, coeffLower, coeffUpper] = ...
%               fit_setup_2exp (direction, varargin)
% Explanation:
%   This function is for fitting equations of the form
%       'a*(1-exp(-x/b))+c*(1-exp(-x/d))' ('rising')
%       'a*exp(-x/b)+c*exp(-x/d)' ('falling')
% Example(s):
%       fit_setup_2exp('rising');
%       fit_setup_2exp('falling', 'MaxScalingFactor', 10);
%       fit_setup_2exp('custom', 'EquationForm', 'a*exp(-x/b)+c*exp(-x/d)+1');
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
% Arguments:    
%       direction   - the direction of the curve
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'rising'  - rising curve
%                       'falling' - falling curve
%                       'custom'  - a custom equation form is provided
%       varargin    - 'AmplitudeEstimate': an estimated amplitude
%                   must be a numeric scalar
%                   default == 1
%                   - 'MaxScalingFactor': maximum scaling factor for amplitude
%                   must be a nonnegative scalar
%                   default == 10
%                   - 'Tau1Init': estimated time constant for the 1st component
%                   must be a nonnegative scalar
%                   default == 10
%                   - 'Tau2Init': estimated time constant for the 2nd component
%                   must be a nonnegative scalar
%                   default == 1
%                   - 'Tau1Range': time constant range for the 1st component
%                   must be a numeric 2-element vector
%                   default == [0, Inf]
%                   - 'Tau2Range': time constant range for the 2nd component
%                   must be a numeric 2-element vector
%                   default == [0, Inf]
%                   - 'EquationForm': equation form if the direction is 'custom'
%                   must be a character vector
%                   default == ''
%
% Used by:    
%       cd/fit_2exp.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2018-10-14 Changed default tau1Init to be 10
% 2018-10-14 Added the possibility of providing a custom equation form
% 

%% Hard-coded parameters
validDirections = {'rising', 'falling', 'custom'};
eqFormRising = 'a*(1-exp(-x/b))+c*(1-exp(-x/d))';
                            % double exponential equation form for rising phase
eqFormFalling = 'a*exp(-x/b)+c*exp(-x/d)';
                            % double exponential equation form for falling phase

%% Default values for optional arguments
amplitudeEstimateDefault = 1;
maxScalingFactorDefault = 10;
tau1InitDefault = 10;
tau2InitDefault = 1;
tau1RangeDefault = [0, Inf];
tau2RangeDefault = [0, Inf];
eqFormDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'direction', ...
    @(x) any(validatestring(x, validDirections)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AmplitudeEstimate', amplitudeEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MaxScalingFactor', maxScalingFactorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Tau1Init', tau1InitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Tau2Init', tau2InitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Tau1Range', tau1RangeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'Tau2Range', tau2RangeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'EquationForm', eqFormDefault, ...
    @(x) validateattributes(x, {'char'}, {'2d'}));
addParameter(iP, 'EquationFormLatex', eqFormDefault, ...
    @(x) validateattributes(x, {'char'}, {'2d'}));

% Read from the Input Parser
parse(iP, direction, varargin{:});
direction = validatestring(direction, validDirections);
amplitudeEstimate = iP.Results.AmplitudeEstimate;
maxScalingFactor = iP.Results.MaxScalingFactor;
tau1Init = iP.Results.Tau1Init;
tau2Init = iP.Results.Tau2Init;
tau1Range = iP.Results.Tau1Range;
tau2Range = iP.Results.Tau2Range;
eqFormCustom = iP.Results.EquationForm;

% Custom equation form will be ignored if direction is not 'custom'
if ~isempty(eqFormCustom) && ~strcmpi(direction, 'custom')
    fprintf(['Warning: Provided equation form will be ignored ', ...
                'unless direction is set to ''custom''!!\n']);
end

% If direction is not 'custom', an equation form must be provided
if strcmpi(direction, 'custom')
    if isempty(eqFormCustom)
        error('Equation form must be provided if direction is ''custom''!');
    elseif ~all(cellfun(@(x) any(strfind(eqFormCustom, x)), ...
                {'a', 'b', 'c', 'd'}))
        error('Equation form must contain ''a'', ''b'', ''c'', ''d''!');
    end
end

%% Preparation
% Decide on equation form
switch direction
case 'rising'
    eqForm = eqFormRising;
case 'falling'
    eqForm = eqFormFalling;
case 'custom'
    eqForm = eqFormCustom;
otherwise
    error('direction unrecognized!!');
end

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
    case {'a' , 'c'}                % amplitudes of each component
        coeffInit(iCoeff) = amplitudeEstimate;
        if amplitudeEstimate > 0
            coeffLower(iCoeff) = 0;
            coeffUpper(iCoeff) = maxScalingFactor * amplitudeEstimate;
        else
            coeffLower(iCoeff) = maxScalingFactor * amplitudeEstimate;
            coeffUpper(iCoeff) = 0;
        end
    case 'b'                        % time constant of first component
        coeffInit(iCoeff) = tau1Init;
        coeffLower(iCoeff) = tau1Range(1);
        coeffUpper(iCoeff) = tau1Range(2);
    case 'd'                        % time constant of second component
        coeffInit(iCoeff) = tau2Init;
        coeffLower(iCoeff) = tau2Range(1);
        coeffUpper(iCoeff) = tau2Range(2);
    otherwise
        error('coeffName is unrecognized!!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if strcmpi(direction, 'rising')
    eqForm = eqFormRising;
elseif strcmpi(direction, 'falling')
    eqForm = eqFormFalling;
else
    error('direction unrecognized!!');
end

eqFormRisingLatex = '$$a(1-e^{-x/b})+c(1-e^{-x/d})$$';
                            % double exponential equation form for rising phase
eqFormFallingLatex = '$$ae^{-x/b}+ce^{-x/d}$$';
                            % double exponential equation form for falling phase

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%