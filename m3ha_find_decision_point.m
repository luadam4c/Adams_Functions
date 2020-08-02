function [indDecision, slopeAtDecision, concavityAtDecision, maxConcavity] = ...
                m3ha_find_decision_point (itm2hDiff, varargin)
%% Finds the indices for the decision points in the 
% Usage: [indDecision, slopeAtDecision, concavityAtDecision, maxConcavity] = ...
%               m3ha_find_decision_point (itm2hDiff, siMs (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       indDecision     - TODO: Description of indDecision
%                       specified as a TODO
%
% Arguments:
%       itm2hDiff   - TODO: Description of itm2hDiff
%                   must be a TODO
%       siMs        - (opt) sampling interval in ms
%                           If not provided, 'tVecs' must be provided
%                   must be a positive vector
%                   default == computed from tVecs
%       varargin    - 'tVecs': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   default == created from siMs and iVecs
%                   - 'FiltWidthMs': filter window width in ms
%                   must be a numeric vector
%                   default == 15 ms
%                   - 'Itm2hDiffLowerLimit': lower limit of discrepancy
%                   must be a numeric vector
%                   default == 1e-9
%                   - 'Itm2hDiffLeftBound': lower limit of discrepancy 
%                                               for LTS region
%                   must be a numeric vector
%                   default == 1e-7
%                   - 'OnlyIfReached': whether to return NaN if decision point 
%                                       is not reached
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/apply_to_all_cells.m
%       cd/argfun.m
%       cd/compute_derivative_trace.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/extract_fields.m
%       cd/extract_subvectors.m
%       cd/find_zeros.m
%       cd/match_time_info.m
%       cd/movingaveragefilter.m
%       cd/parse_peaks.m
%       cd/restrict_values.m
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_simulate_population.m

% File History:
% 2020-05-13 Adapted from code in m3ha_plot_simulated_traces.m
% 2020-05-14 Added 'OnlyIfReached' as an optional argument
% 2020-05-15 Updated algorithm
% 2020-05-16 Changed default filtWidthMs from 30 ms to 15 ms

%% Hard-coded parameters

%% Default values for optional arguments
siMsDefault = [];               % set later
tVecsDefault = [];              % set later
filtWidthMsDefault = 15;        % default filter width is 15 ms
itm2hDiffLowerLimitDefault = 1e-9;
itm2hDiffLeftBoundDefault = 1e-7;
onlyIfReachedDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'itm2hDiff', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['itm2hDiff must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'siMs', siMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'FiltWidthMs', filtWidthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'Itm2hDiffLowerLimit', itm2hDiffLowerLimitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Itm2hDiffLeftBound', itm2hDiffLeftBoundDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'OnlyIfReached', onlyIfReachedDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, itm2hDiff, varargin{:});
siMs = iP.Results.siMs;
tVecs = iP.Results.tVecs;
filtWidthMs = iP.Results.FiltWidthMs;
itm2hDiffLowerLimit = iP.Results.Itm2hDiffLowerLimit;
itm2hDiffLeftBound = iP.Results.Itm2hDiffLeftBound;
onlyIfReached = iP.Results.OnlyIfReached;

%% Preparation
% Count the number of samples for each vector
nSamples = count_samples(itm2hDiff);

% Compute sampling interval(s) and create time vector(s)
if isempty(siMs) && isempty(tVecs)
    error('One of siMs and tVecs must be provided!');
else
    [tVecs, siMs] = ...
        match_time_info(tVecs, siMs, nSamples, 'TimeUnits', 'ms');
end

%% Prepare vectors
% Force itm2hDiff to be above itm2hDiffLowerLimit
itm2hDiff = restrict_values(itm2hDiff, 'LowerBound', itm2hDiffLowerLimit);

% Compute x = the logarithm of m2hDiff
logItm2hDiff = apply_to_all_cells(@log10, itm2hDiff);

% Smooth x over filtWidthMs
logItm2hDiffSmooth = movingaveragefilter(logItm2hDiff, filtWidthMs, siMs);

% Compute dx/dt
%   Note: this is in k(s^-1)
[dxdtVecs, t1Vecs] = compute_derivative_trace(logItm2hDiffSmooth, tVecs);

% Smooth dx/dt over filtWidthMs
dxdtVecs = movingaveragefilter(dxdtVecs, filtWidthMs, siMs);

% Compute d2x/dt2
%   Note: this is in M(s^-2)
[d2xdt2Vecs, t2Vecs] = compute_derivative_trace(dxdtVecs, t1Vecs);

% Smooth d2x/dt2 over filtWidthMs
d2xdt2Vecs = movingaveragefilter(d2xdt2Vecs, filtWidthMs, siMs);

%% Restrict to the rising phase of the maximum peak of itm2hDiff
% Parse maximum peak (considering all peaks) from itm2hDiff, 
%   using itm2hDiffLeftBound as the peak lower bound
peakParams = vecfun(@(x) parse_peaks(x, 'ParseMode', 'maxOfAll', ...
                        'PeakLowerBound', itm2hDiffLeftBound), ...
                    itm2hDiff, 'UniformOutput', true);

% Extract index peak starts and ends
[idxPeakStart, idxPeak] = argfun(@(x) extract_fields(peakParams, x), ...
                                    'idxPeakStart', 'idxPeak');

% Restrict to just the rising phase of the maximum peak
[tRise, t1Rise, t2Rise, logItm2hDiffRise, dxdtRise, d2xdt2Rise] = ...
    argfun(@(x) extract_subvectors(x, 'IndexStart', idxPeakStart, ...
                                    'IndexEnd', idxPeak), ...
            tVecs, t1Vecs, t2Vecs, logItm2hDiff, dxdtVecs, d2xdt2Vecs);

%% Find the decision point
% Find the decision point in each LTS region
[indDecisionRise, slopeAtDecision, concavityAtDecision, maxConcavity] = ...
    m3ha_find_decision_point_helper(logItm2hDiffRise, dxdtRise, ...
                                    d2xdt2Rise, onlyIfReached);

% Convert to the original index
indDecision = idxPeakStart - 1 + indDecisionRise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indDecision, slopeAtDecision, concavityAtDecision, maxConcavity] = ...
                m3ha_find_decision_point_helper (xVecs, dxdtVecs, ...
                                                d2xdt2Vecs, onlyIfReached)

%% Compute index of the maximum value of x 
% Extract the index of maximum value for x 
[~, indMaxValue] = extract_elements(xVecs, 'max');

%% Compute index of the first minimum concavity for x 
%   before the maximum of x

% Restrict to the part of d2xdt2Vecs until the maximum of x is reached
d2xdt2VecsLeft1 = extract_subvectors(d2xdt2Vecs, 'IndexEnd', indMaxValue - 1);

% Extract the index of the first minimum concavity 
%   before the maximum of x is reached
troughParams = vecfun(@(v) parse_peaks(-v, 'ParseMode', 'first'), ...
                        d2xdt2VecsLeft1);
ind2FirstMinConcavityBeforeMax = extract_fields(troughParams, 'idxPeak');
minConcavityBeforeMax = -extract_fields(troughParams, 'peakAmp');

%% Compute index of the maximum concavity of x 
%   between the first minimum concavity of x and the the maximum of x

% Restrict to the part of d2xdt2Vecs after the first minimum
d2xdt2VecsMiddle = extract_subvectors(d2xdt2VecsLeft1, ...
                        'IndexStart', ind2FirstMinConcavityBeforeMax);

% Extract the index of maximum concavity before the maximum of x
%   but before the first minimum concavity
[maxConcavity, ind3MaxConcavityBeforeMax] = ...
    extract_elements(d2xdt2VecsMiddle, 'max');

% Convert to original index
ind2MaxConcavityBeforeMax = ...
    (ind2FirstMinConcavityBeforeMax - 1) + ind3MaxConcavityBeforeMax;

%% Compute index of the last zero concavity of x 
% Compute the index of zero concavity for x 
%   between the first minimum concavity of x 
%       and the maximum concavity of x before the maximum of x

% Restrict to the part before the maximum is reached
d2xdt2VecsRise = extract_subvectors(d2xdt2VecsMiddle, ...
                    'IndexEnd', ind3MaxConcavityBeforeMax);

% Find the last index with d2x/dt2 closest to zero (but positive)
%   between the first minimum concavity of x 
%       and the maximum concavity of x before the maximum of x
ind3LastZeroBeforeMaxConcavityBeforeMax = ...
    find_zeros(d2xdt2VecsRise, 1, 'last', 'IndexChoice', 'positive');

% Convert to original index
ind2LastZeroBeforeMaxConcavityBeforeMax = ...
    (ind2FirstMinConcavityBeforeMax - 1) + ...
        ind3LastZeroBeforeMaxConcavityBeforeMax;

%% Decide on the decision point
% The decision point is either:
%   (1) The last zero before the maximum of d2x/dt2 before the maximum of x
%               or 
%   (2) The first minimum of d2x/dt2 before the maximum of x
%          if this minimum is nonnegative
%               or
%   (3) The maximum of d2x/dt2 before the maximum of x otherwise
%           (not reached)

% If concavity before the maximum reaches negative, 
%   don't use the first minimum as a decision point
ind2FirstMinConcavityBeforeMax(minConcavityBeforeMax < 0) = NaN;

% Use the minimum function 
%   Note: NaN will be ignored unless all are NaN
if onlyIfReached
    ind2Decision = min([ind2LastZeroBeforeMaxConcavityBeforeMax, ...
                        ind2FirstMinConcavityBeforeMax], [], 2);
else
    ind2Decision = min([ind2LastZeroBeforeMaxConcavityBeforeMax, ...
                        ind2FirstMinConcavityBeforeMax, ...
                        ind2MaxConcavityBeforeMax], [], 2);
end

% Find the corresponding indices in x and dx/dt
ind1Decision = ind2Decision;
indDecision = ind2Decision + 1;

% Extract the slope value at the decision point
slopeAtDecision = extract_elements(dxdtVecs, 'specific', 'Index', ind1Decision);

% Extract the concavity value at the decision point
concavityAtDecision = extract_elements(d2xdt2Vecs, 'specific', ...
                                    'Index', ind2Decision);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%