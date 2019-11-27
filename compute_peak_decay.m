function [peakDecaySamples, idxPeakDecay, peakDecayValue] = ...
                compute_peak_decay (vectors, idxPeak, varargin)
%% Computes the peak decays
% Usage: [peakDecaySamples, idxPeakDecay] = ...
%               compute_peak_decay (vectors, idxPeak, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       peakDecaySamples    - peak half widths in samples
%                           specified as a nonnegative vector
%       idxPeakDecay        - indices of start and end of half widths
%                           specified as a positive integer vector
%
% Arguments:
%       vectors     - vectors with peaks
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%       idxPeak     - peak index of each vector
%                   must be a positive integer vector
%       varargin    - 'BaseValue': baseline value of each vector
%                   must be empty or a numeric vector
%                   default == first value of each vector
%                   - 'DecayMethod': method for alignment/truncation
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'exp'  - exponential decay to one time constant
%                       '2exp' - double exponential decay to exp(-1)
%                       '90%'  - decay to 90% of original
%                   default == 'exponential'
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_custom.m
%       cd/fit_2exp.m
%       cd/force_column_vector.m
%       cd/isnumericvector.m
%       cd/match_positions.m
%
% Used by:
%       cd/parse_pulse_response.m

% File History:
% 2018-12-23 Modified from compute_peak_halfwidth.m
% 

%% Hard-coded parameters
validDecayMethods = {'exp', '2exp', '90%'};

%% Default values for optional arguments
baseValueDefault = [];             % set later
decayMethodDefault = '2exp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vectors', ...
    @(x) assert(isnumeric(x) || iscell(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array!']));
addRequired(iP, 'idxPeak', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BaseValue', baseValueDefault, ...
    @(x) assert(isnumericvector(x), ...
                'BaseValue must be either empty or a numeric vector!'));
addParameter(iP, 'DecayMethod', decayMethodDefault, ...
    @(x) any(validatestring(x, validDecayMethods)));

% Read from the Input Parser
parse(iP, vectors, idxPeak, varargin{:});
baseValue = iP.Results.BaseValue;
decayMethod = validatestring(iP.Results.DecayMethod, validDecayMethods);

%% Preparation
% Set default baseline value
%   TODO: Improve this
if isempty(baseValue)
    % Use the first element of each vector
    baseValue = extract_elements(vectors, 'first');
end

% Make sure inputs are in the desired form
[vectors, idxPeak, baseValue] = ...
    argfun(@force_column_vector, vectors, idxPeak, baseValue);

% Count the number of peaks
nPeaks = length(idxPeak);

% Match the counts
[vectors, baseValue] = ...
    argfun(@(x) match_dimensions(x, [nPeaks, 1]), vectors, baseValue);

%% Do the job
% Extract the parts of each vector starting from idxPeak
afterPeak = extract_subvectors(vectors, 'IndexStart', idxPeak);

% Extract the peak value
peakValue = extract_elements(vectors, 'specific', 'Index', idxPeak);

% Decide on the directionFactor
directionFactor = sign(peakValue - baseValue);

% Compute the peak decay in samples
switch decayMethod
    case {'2exp', 'exp'}
        % Compute the peak amplitude
        peakAmplitude = peakValue - baseValue;

        % Compute the expected peak amplitude at one time constant 
        peakDecayAmp = peakAmplitude .* exp(-1);

        % Compute the expected value at one time constant 
        peakDecayValue = baseValue + peakDecayAmp;

        % Shift each afterPeak vector relative to baseline
        % TODO: Make function add_to_vectors.m
        afterPeakShifted = cellfun(@(x, y) x - y, afterPeak, ...
                                num2cell(baseValue), 'UniformOutput', false);

        % Fit each shifted trace after the peak to a double exponential
        %   to extract a decay time constant
        nTraces = numel(afterPeak);
        peakDecaySamples = nan(nTraces, 1);
        parfor iTrace = 1:nTraces
            % Extract things for this trace
            traceToFit = afterPeakShifted{iTrace};
            peakAmplitudeThis = peakAmplitude(iTrace);
            peakDecayAmpThis = peakDecayAmp(iTrace);
            directionFactorThis = directionFactor(iTrace);

            % Only do anything if there is a peak
            if ~isnan(peakAmplitudeThis)
                % Initiate fitObject for parfor
              	fitObject = [];

                % Count the number of samples
                nSamples = numel(traceToFit);

                % Make sure there are enough samples for the fit
                if nSamples > 4
                    % Create an x vector
                    xVec = create_indices([1, nSamples]);

                    % Fit the trace
                    switch decayMethod
                        case '2exp'
                            % Fit this trace to a double exponential
                            [~, fitObject] = ...
                                fit_2exp(traceToFit, 'XVector', xVec, ...
                                        'Direction', 'falling', ...
                                        'AmplitudeEstimate', peakAmplitudeThis);
                        case 'exp'
                            % Fit this trace to a single exponential
                            [~, fitObject] = ...
                                fit_exp(traceToFit, 'XVector', xVec, ...
                                        'Direction', 'falling', ...
                                        'AmplitudeEstimate', peakAmplitudeThis);
                        otherwise
                    end

                    % Evaluate the fit at all points
                    fitVec = fitObject(xVec);

                    % Find the first index that reaches peakDecayAmpThis
                    peakDecaySamples(iTrace) = ...
                        find_custom(fitVec .* directionFactorThis <= ...
                                    peakDecayAmpThis .* directionFactorThis, ...
                                    1, 'first', 'ReturnNaN', true) - 1;
                end
            end
        end
    case '90%'
        % Compute the expected value at 90% decay
        peakDecayValue = baseValue * 0.9 + peakValue * 0.1;

        % Find the first index that reaches the value at 90% decay
        %   This is the 90% decay time in samples
        peakDecaySamples = ...
            cellfun(@(x, y, z) find_custom(x * z <= y * z, 1, 'first', ...
                                                        'ReturnNaN', true), ...
                    afterPeak, num2cell(peakDecayValue), ...
                    num2cell(directionFactor)) - 1;
    otherwise
        error('Code logic error!!')
end

% Compute the index at peak decay in the original vector(s)
idxPeakDecay = idxPeak + peakDecaySamples;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

peakValue = extract_subvectors(vectors, ...
                'EndPoints', transpose([idxPeak, idxPeak]));

[vectors, idxPeak, baseValue] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
            vectors, idxPeak, baseValue);

peakDecaySamples = cellfun(@(x) x(2) - x(1), indHalfWidth);

peakValue = cellfun(@(x, y) x(y), vectors, num2cell(idxPeak));

[idxTemp1, idxTemp2] = ...
    argfun(@(w) cellfun(@(x, y, z) find(x * y >= z * y, 1, 'first'), ...
                    w, num2cell(directionFactor), num2cell(halfPeakValue)), ...
            beforePeakReversed, afterPeak);

% Find the value of that decay
peakDecayValue(iTrace) = 
peakDecayValue = nan(nTraces, 1);

peakShortTimeConstantSamples = nan(nTraces, 1);
peakLongTimeConstantSamples = nan(nTraces, 1);

% Extract coefficient names and values
coeffNames = fitParams.coeffNames;
coeffValues = fitParams.coeffValues;

% Extract the time constants in samples
tau1 = match_positions(coeffValues, coeffNames, 'b'); 
tau2 = match_positions(coeffValues, coeffNames, 'd'); 

% Compute the time constants in samples
peakLongTimeConstantSamples(iTrace) = round(max([tau1, tau2]));
peakShortTimeConstantSamples(iTrace) = round(min([tau1, tau2]));

[fitParams, fitObject] = ...
    fit_2exp(traceToFit, 'Direction', 'falling', ...
            'AmplitudeEstimate', peakAmplitudeThis);

% Extract coefficient names and values
equationStr = fitParams.equationStr;

% Solve for the peak decay in samples (may not be an integer)
peakDecayTemp = ...
    solve_function_at_value(equationStr, peakDecayAmpThis);

% Round to nearest samples
peakDecaySamples(iTrace) = round(peakDecayTemp);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
