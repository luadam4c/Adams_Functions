function [indDecision, slopeValues] = ...
                m3ha_find_decision_point (itm2hDiff, varargin)
%% Finds the indices for the decision points in the 
% Usage: [indDecision, slopeValues] = ...
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
%                   default == 30 ms
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

%% Hard-coded parameters

%% Default values for optional arguments
siMsDefault = [];               % set later
tVecsDefault = [];              % set later
filtWidthMsDefault = 30;            % default filter width is 30 ms
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

%% Prepare itm2hDiff
% Force itm2hDiff to be above itm2hDiffLowerLimit
itm2hDiff = restrict_values(itm2hDiff, 'LowerBound', itm2hDiffLowerLimit);

%% Find the rising phase of the LTS
% Parse maximum peak (considering all peaks) from itm2hDiff, 
%   using itm2hDiffLeftBound as the peak lower bound
peakParams = vecfun(@(x) parse_peaks(x, 'ParseMode', 'maxOfAll', ...
                        'PeakLowerBound', itm2hDiffLeftBound), ...
                    itm2hDiff, 'UniformOutput', true);

% Extract index peak starts and ends
[idxPeakStart, idxPeak] = argfun(@(x) extract_fields(peakParams, x), ...
                                    'idxPeakStart', 'idxPeak');

% Restrict to just the rising phase of the LTS
[tVecsRise, itm2hDiffRise] = ...
    argfun(@(x) extract_subvectors(x, 'IndexStart', idxPeakStart, ...
                                    'IndexEnd', idxPeak), ...
            tVecs, itm2hDiff);

%% Find the decision point
% Find the decision point in each LTS region
[indDecisionRise, slopeValues] = ...
    m3ha_find_decision_point_helper(tVecsRise, itm2hDiffRise, ...
                                    filtWidthMs, siMs, onlyIfReached);

% Convert to the original index
indDecision = idxPeakStart - 1 + indDecisionRise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indDecision, slopeValues] = ...
                m3ha_find_decision_point_helper (tVecs, itm2hDiff, ...
                                            filtWidthMs, siMs, onlyIfReached)

% Compute x = the logarithm of m2hDiff
logItm2hDiff = apply_to_all_cells(@log10, itm2hDiff);

% Smooth x over filtWidthMs
logItm2hDiffSmooth = movingaveragefilter(logItm2hDiff, filtWidthMs, siMs);

% Compute dx/dt
[dxdtVecs, t1Vecs] = compute_derivative_trace(logItm2hDiffSmooth, tVecs);

% Smooth dx/dt over filtWidthMs
dxdtVecs = movingaveragefilter(dxdtVecs, filtWidthMs, siMs);

% Compute d2x/dt2
d2xdt2Vecs = compute_derivative_trace(dxdtVecs, t1Vecs);

% Smooth d2x/dt2 over filtWidthMs
d2xdt2Vecs = movingaveragefilter(d2xdt2Vecs, filtWidthMs, siMs);

% Compute index of maximum concavity for logItm2hDiff 
%   before the maximum of logItm2hDiff
% Extract the index of maximum value for logItm2hDiff 
[~, indMaxValue] = extract_elements(logItm2hDiff, 'max');

% Restrict to the part of logItm2hDiff before the maximum is reached
d2xdt2VecsLeft1 = extract_subvectors(d2xdt2Vecs, ...
                                    'IndexEnd', indMaxValue - 1);

% Extract the index of maximum concavity before the maximum of logItm2hDiff
[~, ind2MaxConcavityBeforeMax] = extract_elements(d2xdt2VecsLeft1, 'max');

% Compute index of zero concavity for logItm2hDiff 
%   before the maximum concavity of logItm2hDiff
% Restrict to the part before the maximum is reached
d2xdt2VecsLeft2 = extract_subvectors(d2xdt2Vecs, ...
                    'IndexEnd', ind2MaxConcavityBeforeMax);

% Find the last index with d2x/dt2 closest to zero 
%   before the maximum is reached
ind2LastZeroBeforeMaxConcavityBeforeMax = ...
    find_zeros(d2xdt2VecsLeft2, 1, 'last');

% The decision point is the last zero before the maximum of d2x/dt2 
%               before maximum of x 
%               or the maximum of d2x/dt2 before maximum of x
%                   if the former doesn't exist
if onlyIfReached
    ind2Decision = ind2LastZeroBeforeMaxConcavityBeforeMax;
else
    ind2Decision = min([ind2MaxConcavityBeforeMax, ...
                        ind2LastZeroBeforeMaxConcavityBeforeMax], [], 2);
end

% Find the corresponding indices in x and dx/dt
ind1Decision = ind2Decision;
indDecision = ind2Decision + 1;

% Extract the slope value at the decision point
slopeValues = extract_elements(dxdtVecs, 'specific', 'Index', ind1Decision);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%