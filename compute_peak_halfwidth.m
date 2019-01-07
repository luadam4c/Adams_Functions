function [halfWidthSamples, indHalfWidthEnds, halfPeakValue] = ...
                compute_peak_halfwidth (vectors, idxPeak, varargin)
%% Computes the half widths for peaks
% Usage: [halfWidthSamples, halfPeakValue, indHalfWidthEnds] = ...
%               compute_peak_halfwidth (vectors, idxPeak, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       halfWidthSamples    - peak half widths in samples
%                           specified as a nonnegative vector
%       halfPeakValue       - value at half peak
%                           specified as a numeric vector
%       indHalfWidthEnds    - indices of start and end of half widths
%                           specified as a cell array of 2-element vectors
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
%       cd/force_column_vector.m
%       cd/isnumericvector.m
%
% Used by:
%       cd/parse_pulse_response.m

% File History:
% 2018-12-22 Created by Adam Lu
% 2018-12-28 Fixed fliplr() -> flipud()
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
addParameter(iP, 'BaseValue', baseValueDefault, ...
    @(x) assert(isnumericvector(x), ...
                'BaseValue must be either empty or a numeric vector!'));

% Read from the Input Parser
parse(iP, vectors, idxPeak, varargin{:});
baseValue = iP.Results.BaseValue;

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
% Extract the parts of each vector ending in idxPeak
beforePeak = extract_subvectors(vectors, 'IndexEnd', idxPeak);

% Reverse beforePeak
beforePeakReversed = cellfun(@flipud, beforePeak, 'UniformOutput', false);

% Extract the parts of each vector starting from idxPeak
afterPeak = extract_subvectors(vectors, 'IndexStart', idxPeak);

% Extract the peak value
peakValue = extract_elements(vectors, 'specific', 'Index', idxPeak);

% Compute the value at half peak
halfPeakValue = mean([baseValue, peakValue], 2);

% Decide on the directionFactor
directionFactor = sign(peakValue - halfPeakValue);

% Find the first index that reaches the value at half peak
[idxTemp1, idxTemp2] = ...
    argfun(@(w) cellfun(@(x, y, z) find_custom(x * z <= y * z, ...
                                            1, 'first', 'ReturnNaN', true), ...
            w, num2cell(halfPeakValue), num2cell(directionFactor)), ...
            beforePeakReversed, afterPeak);

% Compute the endpoints for the half width
idxHalfWidthStart = idxPeak - idxTemp1 + 1;
idxHalfWidthEnd = idxPeak + idxTemp2 - 1;

% Compute the half width in samples
halfWidthSamples = idxHalfWidthEnd - idxHalfWidthStart;

% Output the endpoints for the half width
indHalfWidthEnds = arrayfun(@(x, y) [x; y], idxHalfWidthStart, ...
                            idxHalfWidthEnd, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

peakValue = extract_subvectors(vectors, ...
                'EndPoints', transpose([idxPeak, idxPeak]));

[vectors, idxPeak, baseValue] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
            vectors, idxPeak, baseValue);

halfWidthSamples = cellfun(@(x) x(2) - x(1), indHalfWidth);

peakValue = cellfun(@(x, y) x(y), vectors, num2cell(idxPeak));

[idxTemp1, idxTemp2] = ...
    argfun(@(w) cellfun(@(x, y, z) find(x * y >= z * y, 1, 'first'), ...
                    w, num2cell(directionFactor), num2cell(halfPeakValue)), ...
            beforePeakReversed, afterPeak);

% Compute the end points of each vector ending in idxPeak
endPointsBeforePeak = transpose([ones(nPeaks, 1), idxPeak]);

% Extract the parts of each vector ending in idxPeak
beforePeak = extract_subvectors(vectors, 'EndPoints', endPointsBeforePeak);

% Compute the end points of each vector starting from idxPeak
endPointsAfterPeak = transpose([idxPeak, Inf(nPeaks, 1)]);

% Extract the parts of each vector starting from idxPeak
afterPeak = extract_subvectors(vectors, 'EndPoints', endPointsAfterPeak);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
