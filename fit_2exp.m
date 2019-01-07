function [fitParams, fitObject, goodnessOfFit, algorithmInfo] = ...
                fit_2exp (yVec, varargin)
%% Fits a double exponential curve to data
% Usage: [fitParams, fitObject, goodnessOfFit, algorithmInfo] = ...
%               fit_2exp (yVec, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       fitParams       - fitted parameters returned by parse_fitobject.m,
%                           as well as added fields:
%                               direction
%                       specified as a structure 
%       fitObject       - the cfit object returned by the fit function
%                       specified as a cfit object
%       goodnessOfFit   - the gof structure returned by the fit function
%                       specified as a gof structure 
%       algorithmInfo   - the output structure returned by the fit function
%                           with added fields (see fit_setup_2exp.m):
%                               equationForm
%                               aFittype
%                               coeffInit
%                               coeffLower
%                               coeffUpper
%                               amplitudeEstimate
%                       specified as a gof structure 
% Arguments:
%       yVec        - y vector
%                   must be a numeric vector
%       varargin    - 'XVector': x vector
%                   must be empty or a numeric vector
%                   default == create_indices([1, numel(yVec)])
%                   - 'Direction': the direction of the curve
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'rising'  - rising curve
%                       'falling' - falling curve
%                       'auto'    - automatically detect 
%                                       whether rising or falling
%                       'custom'  - a custom equation form is provided
%                   default == 'auto'
%                   - 'AmplitudeEstimate': an estimated amplitude
%                   must be a numeric scalar
%                   default == yLast - yFirst
%                   - 'MaxScalingFactor': maximum scaling factor for amplitude
%                   must be a nonnegative scalar
%                   default == 10
%                   - 'Tau1Init': estimated time constant for the 1st component
%                   must be a nonnegative scalar
%                   default == mean([xLast, xFirst])
%                   - 'Tau2Init': estimated time constant for the 2nd component
%                   must be a nonnegative scalar
%                   default == mean([xLast, xFirst]) / tauRatio
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
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/fit_setup_2exp.m
%       cd/force_column_vector.m
%       cd/parse_fitobject.m
%
% Used by:
%       cd/compute_peak_decay.m
%       cd/fit_pulse_response.m

% File History
% 2018-10-10 Created by Adam Lu
% 2018-10-11 Made the first argument fitParams
% 2018-10-14 Added the possibility of providing a custom equation form
% 2018-12-24 Made 'XVector' an optional argument
% 

%% Hard-coded parameters
validDirections = {'rising', 'falling', 'custom', 'auto'};
tauRatio = 10;                      % default tau1 / tau2

%% Default values for optional arguments
xVectorDefault = [];                % set later
directionDefault = 'auto';          % set later
amplitudeEstimateDefault = [];      % set later
maxScalingFactorDefault = 10;
tau1InitDefault = [];               % set later
tau2InitDefault = [];               % set later
tau1RangeDefault = [0, Inf];
tau2RangeDefault = [0, Inf];
eqFormDefault = '';

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
addRequired(iP, 'yVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XVector', xVectorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'Direction', directionDefault, ...
    @(x) any(validatestring(x, validDirections)));
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
addParameter(iP, 'EquationForm', eqFormDefault, ...
    @(x) validateattributes(x, {'char'}, {'2d'}));

% Read from the Input Parser
parse(iP, yVec, varargin{:});
xVec = iP.Results.XVector;
direction = validatestring(iP.Results.Direction, validDirections);
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
    elseif ~all(cellfun(@(x) any(strfind(eqFormCustom, x)), {'a', 'b', 'c', 'd'}))
        error('Equation form must contain ''a'', ''b'', ''c'', ''d''!');
    end
end

%% Preparation
% Create an x vector if not provided
if isempty(xVec)
    xVec = create_indices([1, numel(yVec)]);
end

% Check relationships between arguments
if numel(xVec) ~= numel(yVec)
    error('xVec and yVec must have the same number of elements!!');
end

% Force xVec and yVec to be columns
xVec = force_column_vector(xVec);
yVec = force_column_vector(yVec);

% Extract the first and last y value
if isempty(amplitudeEstimate) || strcmpi(direction, 'auto')
    yFirst = yVec(1);
    yLast = yVec(end);
end

% Estimate a direction if not provided
if strcmpi(direction, 'auto')
    if abs(yFirst) <= abs(yLast)
        direction = 'rising';
    else
        direction = 'falling';
    end
end

% Estimate an amplitude if not provided
if isempty(amplitudeEstimate)
    switch direction
    case 'rising'
        amplitudeEstimate = yLast - yFirst;
    case 'falling'
        amplitudeEstimate = yFirst - yLast;
    case 'custom'
        % Find the index of the maximum absolute y value
        [~, idxMaxAbs] = max(abs(yVec));

        % Find the index of the minimum absolute y value
        [~, idxMinAbs] = min(abs(yVec));
        
        % Use the difference between the values corresponding to the
        %   indices above
        amplitudeEstimate = yVec(idxMaxAbs) - yVec(idxMinAbs);
    otherwise
        error('direction unrecognized!!');
    end
end

% Extract the first and last x value
if isempty(tau1Init) || isempty(tau2Init)
    xFirst = xVec(1);
    xLast = xVec(end);
end

% Estimate initial time constants if not provided
if isempty(tau1Init)
    tau1Init = mean([xLast, xFirst]);
end
if isempty(tau2Init)
    tau2Init = mean([xLast, xFirst]) / tauRatio;
end

% Setup fitting type, initial conditions and boundaries
[equationForm, aFittype, coeffInit, coeffLower, coeffUpper] = ...
    fit_setup_2exp(direction, ...
                    'AmplitudeEstimate', amplitudeEstimate, ...
                    'MaxScalingFactor', maxScalingFactor, ...
                    'Tau1Init', tau1Init, ...
                    'Tau2Init', tau2Init, ...
                    'Tau1Range', tau1Range, ...
                    'Tau2Range', tau2Range, ...
                    'EquationForm', eqFormCustom);

%% Do the job
% Fit data to the fitting type
[fitObject, goodnessOfFit, algorithmInfo] = ...
    fit(xVec, yVec, aFittype, 'StartPoint', coeffInit, ...
        'Lower', coeffLower, 'Upper', coeffUpper);

%% Store outputs
% Parse from the fit object
fitParams = parse_fitobject(fitObject);

% Store direction in fitparams
fitParams.direction = direction;

% Store all intermediates in algorithmInfo
algorithmInfo.equationForm = equationForm;
algorithmInfo.aFittype = aFittype;
algorithmInfo.coeffInit = coeffInit;
algorithmInfo.coeffLower = coeffLower;
algorithmInfo.coeffUpper = coeffUpper;
algorithmInfo.amplitudeEstimate = amplitudeEstimate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if ~iscolumn(xVec)
    xVec = xVec(:);
end
if ~iscolumn(yVec)
    yVec = yVec(:);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
