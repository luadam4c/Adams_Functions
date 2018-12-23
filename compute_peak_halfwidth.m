function [halfWidthSamples, halfPeakValue, indHalfWidth] = ...
                compute_peak_halfwidth (vectors, idxPeak, baseValue, varargin)
%% Computes the half widths for peaks
% Usage: [halfWidthSamples, halfPeakValue, indHalfWidth] = ...
%               compute_peak_halfwidth (vectors, idxPeak, baseValue, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       halfWidthSamples    - peak half widths in samples
%                           specified as a nonnegative vector
%       halfPeakValue       - value at half peak
%                           specified as a numeric vector
%       indHalfWidth        - indices of start and end of half peaks
%                           specified as a cell array of 2-element vectors
% Arguments:
%       vectors     - vectors with peaks
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%       idxPeak     - peak index of each vector
%                   must be a positive integer vector
%       varargin    - 'BaseValue': baseline value of each vector
%                   must be empty or a positive integer vector
%                   default == first value of each vector
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
extract_subvectors
force_column_numeric
%       cd/ispositiveintegervector.m
%
% Used by:
%       cd/parse_pulse_response.m

% File History:
% 2018-12-22 Created by Adam Lu
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
    @(x) assert(isempty(x) || ispositiveintegervector(x), ...
                ['idxPeak must be either empty ', ...
                    'or a positive integer vector!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'baseValue', baseValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
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
    @(x) match_dimensions(x, [nPeaks, 1], vectors, baseValue);

%% Do the job
% Extract the peak value
peakValue = extract_subvectors(vectors, 'EndPoints', transpose([idxPeak, idxPeak]));

% Compute the value at half peak
halfPeakValue = mean([baseValue, peakValue], 2);

% Decide on the directionFactor
directionFactor = sign(halfPeakValue - peakValue);

% Extract the parts of each vector before idxPeak
beforePeak = extract_subvectors(vectors, 'EndPoints', transpose([ones(nPeaks, 1), idxPeak]));

% Reverse beforePeak
beforePeakReversed = cellfun(@fliplr, beforePeak, 'UniformOutput', false);

% Extract the parts of each vector after idxPeak
afterPeak = extract_subvectors(vectors, 'EndPoints', transpose([idxPeak, Inf(nPeaks, 1)]));

% Find the first index that reaches the value at half peak

% Find the last index that reaches the value at half peak

% Compute the indices at half width
indHalfWidth = arrayfun(@(x, y, z) [x - y; ], )

% Compute the half width in samples
halfWidthSamples = cellfun(@(x) x(2) - x(1), indHalfWidth);

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

peakValue = arrayfun(@(x, y) x(y), vectors, idxPeak);

[vectors, idxPeak, baseValue] = ...
    argfun(@(x) force_column_numeric(x, 'IgnoreNonVectors', true), ...
            vectors, idxPeak, baseValue);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%