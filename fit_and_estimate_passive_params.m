function [passiveParams, fitResults, fitObject, goodnessOfFit, algorithmInfo] = ...
                fit_and_estimate_passive_params (tVec, yVec, pulseWidth, ...
                                                    pulseAmplitude, varargin)
%% Uses a given pulse width and amplitude to fit and estimate passive parameters from a current pulse response
% Usage: [passiveParams, fitResults, fitObject, goodnessOfFit, algorithmInfo] = ...
%               fit_and_estimate_passive_params (tVec, yVec, pulseWidth, ...
%                                                   pulseAmplitude, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       passiveParams   - parameters returned by estimate_passive_params.m
%                       specified as a structure
%       fitResults      - fitted parameters returned by fit_pulse_response.m
%                       specified as a structure 
%       fitObject       - the cfit object returned by the fit function
%                       specified as a cfit object
%       goodnessOfFit   - the gof structure returned by the fit function
%                       specified as a gof structure 
%       algorithmInfo   - algorithm info returned by fit_pulse_response.m
%                       specified as a scalar structure 
% Arguments:    
%       tVec        - time vector
%                   must be a numeric vector
%       yVec        - response vector
%                   must be a numeric vector
%       pulseWidth  - pulse width in the same units as time vector
%                   must be a nonnegative scalar
%       pulseAmplitude  - pulse amplitude in the same units as response vector
%                   must be a numeric scalar
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
%                   - 'MembraneCapacitance': specific membrane capacitance 
%                                                   (uF/cm^2)
%                   must be a nonnegative scalar
%                   default == 0.88 uF/cm^2
%                   - 'AxialResistivity': axial resistivity (Ohm-cm)
%                   must be a nonnegative scalar
%                   default == 173 Ohm-cm
%                   - 'SeriesResistance': series (pipette) resistance (MOhm)
%                   must be a nonnegative scalar
%                   default == 10 MOhm
%
% Requires:
%       cd/fit_pulse_response.m
%       cd/estimate_passive_params.m
%       cd/merge_structs.m
%
% Used by:    
%       cd/find_passive_params.m
%       cd/plot_cfit_pulse_response.m

% File History:
% 2018-10-11 Moved code from find_passive_params.m
% 2018-10-14 Added the combined phase
% 

%% Hard-coded parameters
validPhaseNames = {'rising', 'falling', 'combined', 'auto'};

%% Default values for optional arguments
phaseNameDefault = 'auto';
amplitudeEstimateDefault = [];
maxScalingFactorDefault = 10;
tau1InitDefault = [];
tau2InitDefault = [];
tau1RangeDefault = [0, Inf];
tau2RangeDefault = [0, Inf];
membraneCapacitanceDefault = 0.88;  % specific membrane capacitance (uF/cm^2)
axialResistivityDefault = 173;      % axial resistivity (Ohm-cm)
seriesResistanceDefault = 10;       % series (pipette) resistance (MOhm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
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
addRequired(iP, 'pulseAmplitude', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PhaseName', phaseNameDefault, ...
    @(x) any(validatestring(x, validPhaseNames)));
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
addParameter(iP, 'MembraneCapacitance', membraneCapacitanceDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'AxialResistivity', axialResistivityDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'SeriesResistance', seriesResistanceDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, tVec, yVec, pulseWidth, pulseAmplitude, varargin{:});
phaseName = validatestring(iP.Results.PhaseName, validPhaseNames);
amplitudeEstimate = iP.Results.AmplitudeEstimate;
maxScalingFactor = iP.Results.MaxScalingFactor;
tau1Init = iP.Results.Tau1Init;
tau2Init = iP.Results.Tau2Init;
tau1Range = iP.Results.Tau1Range;
tau2Range = iP.Results.Tau2Range;
Cm = iP.Results.MembraneCapacitance;
Ra = iP.Results.AxialResistivity;
Rs = iP.Results.SeriesResistance;

% Check relationships between arguments
if length(tVec) ~= length(yVec)
    error('tVec and yVec must have the same length!!');
end

%% Preparation

%% Do the job
if ~isempty(tVec)
    % Fit pooled rising phase data with double exponential
    [fitResults, fitObject, goodnessOfFit, algorithmInfoFit] = ...
        fit_pulse_response (tVec, yVec, pulseWidth, ...
                               'PhaseName', phaseName, ...
                               'AmplitudeEstimate', amplitudeEstimate, ...
                               'MaxScalingFactor', maxScalingFactor, ...
                               'Tau1Init', tau1Init, ...
                               'Tau2Init', tau2Init, ...
                               'Tau1Range', tau1Range, ...
                               'Tau2Range', tau2Range);

    % Estimate passive parameters from long-pulse coefficients
    [passiveParams, algorithmInfoEstimate] = ...
        estimate_passive_params(fitResults, pulseAmplitude, ...
                                'MembraneCapacitance', Cm, ...
                                'AxialResistivity', Ra, ...
                                'SeriesResistance', Rs);

    % Combine the algorithm info structures
    algorithmInfo = merge_structs(algorithmInfoFit, algorithmInfoEstimate);
else
    fitResults = [];
    fitObject = [];
    passiveParams = [];
    goodnessOfFit = [];
    algorithmInfo = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%