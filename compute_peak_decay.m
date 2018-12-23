function [peakDecaySamples, idxPeakDecay] = ...
                compute_peak_decay (vectors, idxPeak, varargin)
%% Computes the peak decays
% Usage: [peakDecaySamples, idxPeakDecay] = ...
%               compute_peak_decay (vectors, idxPeak, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       peakDecaySamples    - peak half widths in samples
%                           specified as a nonnegative vector
%       idxPeakDecay        - indices of start and end of half widths
%                           specified as a positive integer vector
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
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_custom.m
%       cd/force_column_numeric.m
%       cd/isnumericvector.m
%
% Used by:
%       cd/parse_pulse_response.m

% File History:
% 2018-12-23 Modified from compute_peak_halfwidth.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
baseValueDefault = [];             % set later

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
addParameter(iP, 'baseValue', baseValueDefault, ...
    @(x) assert(isnumericvector(x), ...
                'baseValue must be either empty or a numeric vector!'));

% Read from the Input Parser
parse(iP, vectors, idxPeak, varargin{:});
baseValue = iP.Results.baseValue;

%% Preparation
% Set default baseline value
%   TODO: Improve this
if isempty(baseValue)
    % Use the first element of each vector
    baseValue = extract_elements(vectors, 'first');
end

% Make sure inputs are in the desired form
[vectors, idxPeak, baseValue] = ...
    argfun(@force_column_numeric, vectors, idxPeak, baseValue);

% Count the number of peaks
nPeaks = length(idxPeak);

% Match the counts
[vectors, baseValue] = ...
    argfun(@(x) match_dimensions(x, [nPeaks, 1]), vectors, baseValue);

%% Do the job
% Extract the peak value
peakValue = extract_elements(vectors, 'specific', 'Index', idxPeak);

% Compute the value at half peak
halfPeakValue = mean([baseValue, peakValue], 2);

% Decide on the directionFactor
directionFactor = sign(peakValue - halfPeakValue);

% Compute the end points of each vector before idxPeak
endPointsBeforePeak = transpose([ones(nPeaks, 1), idxPeak]);

% Extract the parts of each vector before idxPeak
beforePeak = extract_subvectors(vectors, 'EndPoints', endPointsBeforePeak);

% Reverse beforePeak
beforePeakReversed = cellfun(@fliplr, beforePeak, 'UniformOutput', false);

% Compute the end points of each vector after idxPeak
endPointsAfterPeak = transpose([idxPeak, Inf(nPeaks, 1)]);

% Extract the parts of each vector after idxPeak
afterPeak = extract_subvectors(vectors, 'EndPoints', endPointsAfterPeak);

% Find the first index that reaches the value at half peak
[idxTemp1, idxTemp2] = ...
    argfun(@(w) cellfun(@(x, y, z) find_custom(x * y >= z * y, 1, 'first', ...
                                                'ReturnNaN', true), ...
            w, num2cell(directionFactor), num2cell(halfPeakValue)), ...
            beforePeakReversed, afterPeak);

% Compute the endpoints for the half width
idxHalfWidthStart = idxPeak - idxTemp1;
idxHalfWidthEnd = idxPeak + idxTemp2;

% Compute the half width in samples
peakDecaySamples = idxHalfWidthEnd - idxHalfWidthStart;

% Output the endpoints for the half width
idxPeakDecay = ...
    arrayfun(@(x, y) [x; y], idxHalfWidthStart, idxHalfWidthEnd, ...
                'UniformOutput', false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

peakValue = extract_subvectors(vectors, ...
                'EndPoints', transpose([idxPeak, idxPeak]));

[vectors, idxPeak, baseValue] = ...
    argfun(@(x) force_column_numeric(x, 'IgnoreNonVectors', true), ...
            vectors, idxPeak, baseValue);

peakDecaySamples = cellfun(@(x) x(2) - x(1), indHalfWidth);

peakValue = cellfun(@(x, y) x(y), vectors, num2cell(idxPeak));

[idxTemp1, idxTemp2] = ...
    argfun(@(w) cellfun(@(x, y, z) find(x * y >= z * y, 1, 'first'), ...
                    w, num2cell(directionFactor), num2cell(halfPeakValue)), ...
            beforePeakReversed, afterPeak);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%