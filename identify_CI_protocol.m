function [isCI, endPointsPulse] = identify_CI_protocol (iVecs, siMs)
%% Identifies whether a set of current vectors is a current injection protocol, and if so, what the range of the current injection is
% Usage: [isCI, endPointsPulse] = identify_CI_protocol (iVecs, siMs)
% Explanation:
%       TODO: Explain strategy
%
% Outputs: TODO
%       isCI            - boolean whether or not abfdata is a current injection
%       endPointsPulse  - values representing the start and end of the injection
%
% Arguments: TODO
%       iVecs       - current vectors
%       siMs        - sampling interval in milliseconds
%
% Requires:
%       cd/freqfilter.m
%       cd/medianfilter.m
%       cd/parse_abf.m TODO
%
% Used by:
%       cd/parse_abf.m
%       cd/plot_FI.m

% File History:
% 2017-04-26 Created by BT
% 2018-07-24 AL - Reorganized comments
% 2018-07-24 AL - Updated usage of identify_channels.m
% 2018-08-02 BT - Removed second derivative in CI range finding and 
%                   changed CI ID parameters
% 2018-08-04 BT - Restored 2nd derivative, added more CI ID parameters
% 2018-09-25 AL - No longer uses identify_channels.m
% 2019-11-05 AL - Changed the 2nd argument from siUs to siMs
% TODO: Make the first argument iVecsORfileName 
%           and use parse_abf.m as in identify_eLFP.m
% TODO: Make siMs an optional argument with default 0.1
% TODO: Add input parser
% TODO: Doesn't work in some cases?
%       See /media/ashleyX/Recordings/20180927/2018_09_27_0001_traces

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
cutoffFreq = 50;            % cutoff frequency (Hz) for lowpass filter
nPoles = 8;                 % poles for lowpass filter
swpSpacingTolerance = 2;    % tolerance for differences in spacing between the sweeps
zeroTolerance = 1;          % tolerance for derivative over current injection range
rangeTolerance = 3;         % tolerance for differences in CI start/end points via 1st/2nd derivs
minCIRangeTolerance = 500;  % minimum range a CI must take

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of sweeps
nSweeps = size(iVecs, 2);

% If there are no sweeps, this is not a current injection protocol
if nSweeps == 0 || isempty(iVecs)
    isCI = false;
    endPointsPulse = NaN;
    return;
end

%% Do the job
% Apply lowpass filter
iVecsLowPass = freqfilter(iVecs, cutoffFreq, siMs / MS_PER_S, ...
                            'FilterType', 'low', 'FilterOrder', nPoles);

% Apply a 3rd order median filter to take out some noisy spikes
iVecsFilt = medianfilter(iVecsLowPass);

% Greatest amplitude CI sweep
[~, maxSwp] = max(mean(abs(iVecsFilt), 1));

% First derivative of current data for the max sweep 
%   (if CI, all CI should be on same time interval)
deriv1Curr = diff(iVecsFilt(:, maxSwp), 1);

% second derivative of current data for the max sweep 
%   (if CI, all CI should be on same time interval)
deriv2Curr = diff(iVecsFilt(:, maxSwp), 2);

% function for second largest val
max2 = @(x) max(x(x < max(x)));

% function for second smallest val
min2 = @(x) min(x(x > min(x)));

% Store potential CI range options
ciOptions = [];

% 8 options for CI range, max & second max, min & second min, 
%   for both 1st & 2nd deriv
% CI has jump up and jump down, using derivatives to find where the jumps are
% Sometimes first derivatives don't catch the jumps
[~, ci1] = max(deriv1Curr);
[~, ci2] = max2(deriv1Curr);
[~, ci3] = min(deriv1Curr);
[~, ci4] = min2(deriv1Curr);
[~, ci5] = max(deriv2Curr);
[~, ci6] = max2(deriv2Curr);
[~, ci7] = min(deriv2Curr);
[~, ci8] = min2(deriv2Curr);
ciOptions = sort([ci1, ci2, ci3, ci4, ci5, ci6, ci7, ci8]);

% Removes possible options near each other
ciOptions = uniquetol(ciOptions, rangeTolerance, 'DataScale', 10);

% AL - Added for a bug
if ciOptions < 2
    isCI = false;
    endPointsPulse = NaN;
    return
end

% With sorted ciOptions, permutations of possible CI ranges
possibleRangeCIs = nchoosek(ciOptions, 2);

% Temporary comparison for maximum average of each permutation
maxRangeAvg = 0;

% Index of possibleRangeCIs for most likely set
likelyRangeSet = 1;

% Each possible CI range
for x = 1:size(possibleRangeCIs, 1)
    % Averages each sweep within set
    swpAvgsInRange = ...
        mean(iVecsFilt(possibleRangeCIs(x, 1):possibleRangeCIs(x, 2), :));

    % Correct injection interval should have greatest average
    if max(abs(swpAvgsInRange)) > abs(maxRangeAvg)
        maxRangeAvg = max(abs(swpAvgsInRange));
        likelyRangeSet = x;
    end
end

% Store most likely set as CI range
idxPulseStart = possibleRangeCIs(likelyRangeSet, 1);
idxPulseEnd = possibleRangeCIs(likelyRangeSet, 2);

% Return pulse endpoints
endPointsPulse = [idxPulseStart, idxPulseEnd];

% Derivatives of all current data over injection range
injectionDataDiff = diff(iVecsFilt(idxPulseStart:idxPulseEnd, :), 1);

% All current data over injection range
injectionDataVals = iVecsFilt(idxPulseStart:idxPulseEnd, :);

% All current data outside injection range
noInjectionDataVals = [iVecsFilt(1:idxPulseStart, :); ...
                        iVecsFilt(idxPulseEnd:end, :)];

% Get averages of each derivative
injectionDataDiffMeans = mean(injectionDataDiff, 1);

% Means of current data over injection range
injectionDataValsMeans = mean(injectionDataVals, 1);

% Means of current data outside injection range
noInjectionDataValsMeans = mean(noInjectionDataVals, 1);

% Difference between sweep derivatives should be minimal
reduction = abs(diff(injectionDataDiffMeans,2));

% Injections should be non-zero, one is tolerated
diffInjection = abs(injectionDataValsMeans - noInjectionDataValsMeans);
diffInjection(find(diffInjection < zeroTolerance, 1)) = [];

% Determine whether it's a current injection protocol
% Criteria:  
%   (1) Are the injection derivatives roughly the same?
%   (2) Are the injection derivatives roughly zero?
%   (3) Are the injections different from the rest of the recording?
%   (4) Is the injection range appropriately sized?
%   (5) Are there multiple sweeps being done?
if reduction < swpSpacingTolerance && ...
        injectionDataDiffMeans < zeroTolerance && ...
        diffInjection > zeroTolerance && ...
        idxPulseEnd - idxPulseStart > minCIRangeTolerance && ...
        nSweeps > 1
    isCI = true;
else
    isCI = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%