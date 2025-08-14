function vector = interpolate_nan_values(vector)
%% Interpolate NaN values in input vector to nearest values
% Usage: vector = interpolate_nan_values(vector)
% Explanation:
%       Replaces NaN values in a vector with interpolated values
%
% Example(s):
%       vector = interpolate_nan_values([3, 5, NaN, NaN, 8, 10])
%       vector = interpolate_nan_values([3, 5, NaN, NaN, 8, 10]')
%
% Outputs:
%       vector  - vector with no NaN values
%               specified as a numeric vector
%
% Arguments:    
%       vector  - vector with or without NaN values
%               specified as a numeric vector
%       
% Requires:
%
% Used by:
%       \Shared\scAAV\analyze_reachr_motion.m

% File History:
% 2025-08-13 Renamed from util_interp1nans.m

% Find the indices with and without NaN
indNan = find(isnan(vector));
indNonNaN = find(~isnan(vector));

% Interpolate at the indices with NaN
valuesInterpolated = interp1(indNonNaN, vector(indNonNaN), indNan);

% Replace NaN with interpolated values
vector(indNan) = valuesInterpolated;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%