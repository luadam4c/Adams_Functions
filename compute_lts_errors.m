function errorStruct = compute_lts_errors (ltsTableSim, ltsTableRec, varargin)
%% Computes low-threshold spike errors for single neuron data
% Usage: errorStruct = compute_lts_errors (ltsTableSim, ltsTableRec, varargin)
% Explanation:
%       Let b = (1/(1+match2FeatureErrorRatio))
%           c = featureWeights(1)/sum(featureWeights)
%           d = featureWeights(2)/sum(featureWeights)
%           e = featureWeights(3)/sum(featureWeights)
%
%       Then:
%       avgLtsError = 
%           ltsMatchError * (1 - b) + ...
%           avgLtsFeatureError * b
%       ltsMatchError = missedLtsError .* sum(hasMissedLts) + ...
%                           falseLtsError .* sum(hasFalseLts);
%       avgLtsFeatureError = ...
%           avgLtsAmpError * c + ...
%           avgLtsDelayError * d + ...
%           avgLtsSlopeError * e
%
%       In more detail:
%       avgLtsError = 
%           ltsMatchError* (1 - b) + ...
%           avgLtsAmpError* b * c + ...
%           avgLtsDelayError* b * d + ...
%           avgLtsSlopeError* b * e
%
% Example(s):
%       TODO
%
% Outputs:
%       errorStruct - a structure of all the errors computed, with fields:
%                       normalizeError (only if normalizeError is true)
%                       normAvgLtsError (only if normalizeError is true)
%                       initLtsError (only if normalizeError is true)
%                       avgLtsError
%                       ltsMatchError
%                       avgLtsFeatureError
%                       avgLtsAmpError
%                       avgLtsDelayError
%                       avgLtsSlopeError
%                       featureWeights
%                       missedLtsError
%                       falseLtsError
%                       match2FeatureErrorRatio
%                       ltsAmpErrors
%                       ltsDelayErrors
%                       ltsSlopeErrors
%                       ltsSweepWeights
%                       hasLtsInBoth
%                       ltsAmpUncertainty
%                       ltsDelayUncertainty
%                       peakPromNormError
%                       peakWidthNormError
%                       slopeUncertainty
%                       baseNoise
%                       sweepWeights
%                   specified as a scalar structure
% Arguments:    
%       ltsTableSim - a table of lts features from simulated voltage traces
%                   must be a table with columns:
%                       ltsPeakTime
%                       ltsPeakValue
%                       maxSlopeValue
%       ltsTableRec - a table of lts features from recorded voltage traces
%                   must be a table with columns:
%                       ltsPeakTime
%                       ltsPeakValue
%                       maxSlopeValue
%       varargin    - 'BaseNoise': baseline noise value(s)
%                   must be a numeric vector
%                   default == none provided
%                   - 'SweepWeights': sweep weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_weighted_average.m
%                   - 'FeatureWeights': LTS feature weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == [2; 2; 2]
%                   - 'MissedLtsError': a dimensionless error that penalizes 
%                                       a misprediction of the existence of LTS
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 2
%                   - 'FalseLtsError': a dimensionless error that penalizes 
%                                       a misprediction of the absence of LTS
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 0.5
%                   - 'Match2FeatureErrorRatio': ratio of LTS match error to 
%                                                   LTS feature error
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 1
%                   - 'NormalizeError': whether to normalize errors 
%                                       by an initial error
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'InitLtsError': initial low-threshold spike errors
%                   must be empty or a numeric vector with length == nSweeps
%                   default == []
%
% Requires:
%       cd/argfun.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/compute_rms_error.m
%       cd/compute_weighted_average.m
%       cd/create_time_vectors.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/iscellnumericvector.m
%       cd/isnumericvector.m
%       cd/match_row_count.m
%       cd/normalize_by_initial_value.m
%
% Used by:
%       cd/compute_single_neuron_errors.m

% File History:
% 2019-11-15 Adapted from compute_sweep_errors.m
% 2019-11-16 Added 'FeatureWeights' as an optional argument
% 2019-11-16 Added 'LtsExistError' as an optional argument
% 2019-11-18 Fixed bug with hasLtsInBoth 
% 2019-11-21 For singleneuronfitting61, fixed bug for the computation of 
%               normalizedDifference. It was fitting with only ltsMatchError
%               in singleneuronfitting60
% 2019-11-25 LTS amplitude uncertainty is now half of peak prominence
% 2019-11-28 Now computes LTS amplitude error assymetrically so that
%               negative errors are penalized 3 times as much
% 2019-11-29 Added ltsMatchError
% 2019-11-29 Now computes LTS amplitude error assymetrically so that
%               positive errors are penalized 10 times less
% 2019-11-29 LTS amplitude uncertainty is now 1/4 of peak prominence
% 2019-12-18 Added 'Match2FeatureErrorRatio' as an optional parameter

%% Hard-coded parameters
% Consistent with singleneuronfitting71.m
defaultLtsFeatureWeights = [2; 2; 2];   % default weights for optimizing 
                                        %   LTS statistics
defaultMissedLtsError = 2;              % how much error (dimensionless) to 
                                        %   penalize a sweep that mispredicted 
                                        %   the existence of an LTS
defaultFalseLtsError = 0.5;             % how much error (dimensionless) to 
                                        %   penalize a sweep that mispredicted 
                                        %   the absence of an LTS
defaultMatch2FeatureErrorRatio = 1;     % default error ratio between
                                        %   match error and avg feature error

%% Default values for optional arguments
baseNoiseDefault = [];          % set later
sweepWeightsDefault = [];       % set in compute_weighted_average.m
featureWeightsDefault = [];     % set later
missedLtsErrorDefault = [];     % set later
falseLtsErrorDefault = [];      % set later
match2FeatureErrorRatioDefault = []; % set later
normalizeErrorDefault = false;  % don't normalize errors by default
initLtsErrorDefault = [];   % no initial error values by default

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
addRequired(iP, 'ltsTableSim', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addRequired(iP, 'ltsTableRec', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BaseNoise', baseNoiseDefault, ...
    @(x) assert(isnumericvector(x), 'BaseNoise must be a numeric vector!'));
addParameter(iP, 'SweepWeights', sweepWeightsDefault, ...
    @(x) assert(isnumericvector(x), 'SweepWeights must be a numeric vector!'));
addParameter(iP, 'FeatureWeights', featureWeightsDefault, ...
    @(x) assert(isnumericvector(x), 'FeatureWeights must be a numeric vector!'));
addParameter(iP, 'MissedLtsError', missedLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'MissedLtsError must be a numeric vector!'));
addParameter(iP, 'FalseLtsError', falseLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'FalseLtsError must be a numeric vector!'));
addParameter(iP, 'Match2FeatureErrorRatio', match2FeatureErrorRatioDefault, ...
    @(x) assert(isnumericvector(x), 'Match2FeatureErrorRatio must be a numeric vector!'));
addParameter(iP, 'NormalizeError', normalizeErrorDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'InitLtsError', initLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'InitLtsError must be a numeric vector!'));

% Read from the Input Parser
parse(iP, ltsTableSim, ltsTableRec, varargin{:});
baseNoise = iP.Results.BaseNoise;
sweepWeights = iP.Results.SweepWeights;
featureWeights = iP.Results.FeatureWeights;
missedLtsError = iP.Results.MissedLtsError;
falseLtsError = iP.Results.FalseLtsError;
match2FeatureErrorRatio = iP.Results.Match2FeatureErrorRatio;
normalizeError = iP.Results.NormalizeError;
initLtsError = iP.Results.InitLtsError;

%% Preparation
% Decide on LTS feature weights
if isempty(featureWeights)
    featureWeights = defaultLtsFeatureWeights;
end

% Decide on missed LTS error
if isempty(missedLtsError)
    missedLtsError = defaultMissedLtsError;
end

% Decide on false LTS error
if isempty(falseLtsError)
    falseLtsError = defaultFalseLtsError;
end

% Decide on LTS match to feature error ratio
if isempty(match2FeatureErrorRatio)
    match2FeatureErrorRatio = defaultMatch2FeatureErrorRatio;
end

% Extract low-threshold spike feature values to compute error for 
%   from both conditions
ltsPeakTimeSim = ltsTableSim.ltsPeakTime;
ltsPeakValueSim = ltsTableSim.ltsPeakValue;
maxSlopeValueSim = ltsTableSim.maxSlopeValue;

ltsPeakTimeRec = ltsTableRec.ltsPeakTime;
ltsPeakValueRec = ltsTableRec.ltsPeakValue;
maxSlopeValueRec = ltsTableRec.maxSlopeValue;

% Extract other feature values for uncertainty calculations
peakWidthSim = ltsTableSim.peakWidth;
peakWidthRec = ltsTableRec.peakWidth;
peakPromSim = ltsTableSim.peakProm;
peakPromRec = ltsTableRec.peakProm;
maxNoiseRec = ltsTableRec.maxNoise;

% Decide on baseline noise
if isempty(baseNoise)
    baseNoise = maxNoiseRec;
end

% Compute default sweep weights
if isempty(sweepWeights)
    % Compute sweep weights based on baseNoise
    sweepWeights = 1 ./ baseNoise;

    % Normalize so that it sums to one
    sweepWeights = sweepWeights / sum(sweepWeights);
end

% Count the number of sweeps
nSweeps = max(height(ltsTableSim), height(ltsTableRec));

% Match row counts for sweep-dependent variables with the number of sweeps
[ltsPeakTimeSim, ltsPeakValueSim, maxSlopeValueSim, ...
        ltsPeakTimeRec, ltsPeakValueRec, maxSlopeValueRec, ...
        baseNoise, sweepWeights] = ...
    argfun(@(x) match_row_count(x, nSweeps), ...
            ltsPeakTimeSim, ltsPeakValueSim, maxSlopeValueSim, ...
            ltsPeakTimeRec, ltsPeakValueRec, maxSlopeValueRec, ...
            baseNoise, sweepWeights);

%% Modify sweep weights for LTS errors
% Initialize LTS sweep weights with sweep weights
ltsSweepWeights = sweepWeights;

% Put the features side by side
ltsFeaturesSim = [ltsPeakTimeSim, ltsPeakValueSim, maxSlopeValueSim];
ltsFeaturesRec = [ltsPeakTimeRec, ltsPeakValueRec, maxSlopeValueRec];

% Determine whether each sweep has an LTS in both simulated and recorded traces
hasLtsInBoth = all(~isnan(ltsFeaturesSim), 2) & all(~isnan(ltsFeaturesRec), 2);

% Make the sweeps without both no LTS weight zero
%   Note: One cannot compute a feature error if either the simulated or 
%           the recorded trace doesn't have an LTS
ltsSweepWeights(~hasLtsInBoth) = 0;

%% Compute the LTS mismatch error
% Determine whether each sweep has an LTS in recorded but not simulated traces
hasMissedLts = any(isnan(ltsFeaturesSim), 2) & all(~isnan(ltsFeaturesRec), 2);

% Determine whether each sweep has an LTS in simulated but not recorded traces
hasFalseLts = all(~isnan(ltsFeaturesSim), 2) & any(isnan(ltsFeaturesRec), 2);

% Compute the LTS mismatch error
ltsMatchError = missedLtsError .* sum(hasMissedLts) + ...
                    falseLtsError .* sum(hasFalseLts);

%% Compute feature uncertainties
% The amplitude uncertainty should be close to 1/4 of peak prominence
ltsAmpUncertainty = abs(peakPromRec ./ 4);

% The peak delay uncertainty should be close to peak width
ltsDelayUncertainty = abs(peakWidthRec);

% Compute the normalized error for peak prominence
peakPromNormError = abs((peakPromSim - peakPromRec) ./ peakPromRec);

% Compute the normalized error for peak width
peakWidthNormError = abs((peakWidthSim - peakWidthRec) ./ peakWidthRec);

% The slope uncertainty can be computed from a propagation of errors
slopeUncertainty = sqrt(peakPromNormError .^ 2 + peakWidthNormError .^ 2) ...
                    .* abs(maxSlopeValueRec);

%% Compute errors
% Compute dimensionless LTS errors for each sweep
ltsAmpErrors = compute_feature_error(ltsPeakValueRec, ltsPeakValueSim, ...
                                            ltsAmpUncertainty);
ltsDelayErrors = compute_feature_error(ltsPeakTimeRec, ltsPeakTimeSim, ...
                                            ltsDelayUncertainty);
ltsSlopeErrors = compute_feature_error(maxSlopeValueRec, maxSlopeValueSim, ...
                                            slopeUncertainty);

% Modify LTS amplitude errors so that positive errors are penalized 
%   10 times less
ltsAmpErrors(ltsAmpErrors >= 0) = ltsAmpErrors(ltsAmpErrors >= 0) ./ 10;

% Compute weighted-root-mean-squared-averaged LTS errors (dimensionless)
[avgLtsAmpError, avgLtsDelayError, avgLtsSlopeError] = ...
    argfun(@(x) compute_weighted_average(x, 'Weights', ltsSweepWeights, ...
                    'IgnoreNaN', true, 'AverageMethod', 'root-mean-square'), ...
            ltsAmpErrors, ltsDelayErrors, ltsSlopeErrors);

% Put the feature errors together
featureErrors = [avgLtsAmpError; avgLtsDelayError; avgLtsSlopeError];

% Compute the weighted average of 
%   the LTS feature errors, weighted by featureWeights
avgLtsFeatureError = compute_weighted_average(featureErrors, ...
                        'Weights', featureWeights, ...
                        'AverageMethod', 'linear', 'IgnoreNaN', true);

% Combine the errors and weights
errorsToAverage = [avgLtsFeatureError; ltsMatchError];
weightsForErrors = [1; match2FeatureErrorRatio];

% Average LTS error (dimensionless) is the weighted average of 
%   the average LTS feature error and the LTS match error, 
%   weighted by match2FeatureErrorRatio
avgLtsError = compute_weighted_average(errorsToAverage, 'IgnoreNan', true, ...
                        'Weights', weightsForErrors, 'AverageMethod', 'linear');

% If requested, make errors dimensionless by 
%   storing or dividing by an initial error value
if normalizeError
    [normAvgLtsError, initLtsError] = ...
        normalize_by_initial_value(avgLtsError, initLtsError);
end

%% Store in output errors structure
% Store normalized errors
if normalizeError
    errorStruct.normalizeError = normalizeError;
    errorStruct.normAvgLtsError = normAvgLtsError;
    errorStruct.initLtsError = initLtsError;
end

% Store other errors in descending order of importance
errorStruct.avgLtsError = avgLtsError;
errorStruct.ltsMatchError = ltsMatchError;
errorStruct.avgLtsFeatureError = avgLtsFeatureError;
errorStruct.avgLtsAmpError = avgLtsAmpError;
errorStruct.avgLtsDelayError = avgLtsDelayError;
errorStruct.avgLtsSlopeError = avgLtsSlopeError;
errorStruct.ltsFeatureWeights = featureWeights;
errorStruct.missedLtsError = missedLtsError;
errorStruct.falseLtsError = falseLtsError;
errorStruct.match2FeatureErrorRatio = match2FeatureErrorRatio;
errorStruct.ltsAmpErrors = ltsAmpErrors;
errorStruct.ltsDelayErrors = ltsDelayErrors;
errorStruct.ltsSlopeErrors = ltsSlopeErrors;
errorStruct.ltsSweepWeights = ltsSweepWeights;
errorStruct.hasLtsInBoth = hasLtsInBoth;
errorStruct.ltsAmpUncertainty = ltsAmpUncertainty;
errorStruct.ltsDelayUncertainty = ltsDelayUncertainty;
errorStruct.peakPromNormError = peakPromNormError;
errorStruct.peakWidthNormError = peakWidthNormError;
errorStruct.slopeUncertainty = slopeUncertainty;
errorStruct.baseNoise = baseNoise;
errorStruct.sweepWeights = sweepWeights;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function featureError = compute_feature_error (valueRec, valueSim, uncertainty)
%% Computes a dimensionless feature error relative to the uncertainty

% Initialize errors
featureError = nan(size(valueRec));

% Put values side by side
allValues = [valueRec, valueSim];

% Determine whether each sweep each condition has no features
noFeature = isnan(allValues);

% Determine whether feature does not exist in either condition
noFeatureInEither = any(noFeature, 2);

% The feature error is not available in this case
featureError(noFeatureInEither) = NaN;

% Determine whether feature exists in both conditions
hasFeatureInBoth = ~noFeatureInEither;

% Compute the normalized difference between simulated and 
%   recorded feature values
normalizedDifference = diff(allValues, 1, 2) ./ abs(uncertainty);

% The feature error is the normalized difference in this case
featureError(hasFeatureInBoth) = normalizedDifference(hasFeatureInBoth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
