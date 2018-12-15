function [isCI, rangeCI] = identify_CI_protocol (iVecs, siUs)
%% Identifies whether a set of current vectors is a current injection protocol, and if so, what the range of the current injection is
% Usage: [isCI, rangeCI] = identify_CI_protocol (iVecs, siUs)
% Explanation:
%       TODO: Explain strategy
%
% Outputs:
%       isCI        - boolean whether or not abfdata is a current injection
%       rangeCI     - values representing the start and end of the injection
% Arguments:
%       iVecs       - current vectors
%       siUs        - sampling interval in microseconds
%
% Requires:
%       cd/parse_abf.m
%
% Used by:
%       cd/parse_abf.m
%       cd/plot_FI.m
%
% File History:
% 2017-04-26 Created by BT
% 2018-07-24 AL - Reorganized comments
% 2018-07-24 AL - Updated usage of identify_channels.m
% 2018-08-02 BT - Removed second derivative in CI range finding and 
%                   changed CI ID parameters
% 2018-08-04 BT - Restored 2nd derivative, added more CI ID parameters
% 2018-09-25 AL - No longer uses identify_channels.m
% TODO: Make the first argument iVecsORfileName 
%           and use parse_abf.m as in identify_eLFP.m
% TODO: Doesn't work in some cases?
%       See /media/ashleyX/Recordings/20180927/2018_09_27_0001_traces

%% Hard-coded parameters
nPoles = 8;                 % poles for lowpass filter
swpSpacingTolerance = 2;    % tolerance for differences in spacing between the sweeps
zeroTolerance = 1;          % tolerance for derivative over current injection range
rangeTolerance = 3;         % tolerance for differences in CI start/end points via 1st/2nd derivs
minCIRangeTolerance = 500;  % minimum range a CI must take

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of sweeps
nSweeps = size(iVecs, 2);

% If there are no sweeps, this is not a current injection protocol
if nSweeps == 0 || isempty(iVecs)
    isCI = false;
    rangeCI = NaN;
    return;
end

% Apply lowpass filter
currentVec = freqfilter(iVecs, 50, siUs * 10^-6, 'FilterType', 'low');

% Apply a 3rd order median filter to take out some noisy spikes
for x = 1:size(currentVec,2)
    currentVec(:, x) = medfilt1(currentVec(:, x));
end

% Greatest amplitude CI sweep
[~, max_swp] = max(mean(abs(currentVec),1));

% First derivative of current data for the max sweep (if CI, all CI should be on same time interval)
deriv1_cd = diff(currentVec(:,max_swp),1);

% second derivative of current data for the max sweep (if CI, all CI should be on same time interval)
deriv2_cd = diff(currentVec(:,max_swp),2);

% function for second largest val
max2 = @(x) max(x(x<max(x)));

% function for second smallest val
min2 = @(x) min(x(x>min(x)));

% Store potential CI range options
ciOptions = [];

% 8 options for CI range, max & second max, min & second min, for both 1st & 2nd deriv
% CI has jump up and jump down, using derivatives to find where the jumps are
% Sometimes first derivatives don't catch the jumps
[~, ci1] = max(deriv1_cd);

[~, ci2] = max2(deriv1_cd);

[~, ci3] = min(deriv1_cd);

[~, ci4] = min2(deriv1_cd);

[~, ci5] = max(deriv2_cd);

[~, ci6] = max2(deriv2_cd);

[~, ci7] = min(deriv2_cd);

[~, ci8] = min2(deriv2_cd);

ciOptions = sort([ci1 ci2 ci3 ci4 ci5 ci6 ci7 ci8]);

% Removes possible options near each other
ciOptions = uniquetol(ciOptions, rangeTolerance, 'DataScale', 10);

% AL - Added for a bug
if ciOptions < 2
    isCI = false;
    rangeCI = NaN;
    return
end

% With sorted ciOptions, permutations of possible CI ranges
possibleRangeCIs = nchoosek(ciOptions, 2);

% Temporary comparison for maximum average of each permutation
maxRangeAvg = 0;

% Index of possibleRangeCIs for most likely set
likelyRangeSet = 1;

% Each possible CI range
for x = 1:size(possibleRangeCIs,1)
    % Averages each sweep within set
    swp_avgs_inrange = mean(currentVec(possibleRangeCIs(x, 1):possibleRangeCIs(x, 2), :));

    % Correct injection interval should have greatest average
    if max(abs(swp_avgs_inrange)) > abs(maxRangeAvg)
        maxRangeAvg = max(abs(swp_avgs_inrange));
        likelyRangeSet = x;
    end
end

% Store likely set as CI range
rangeCI = [possibleRangeCIs(likelyRangeSet,1) possibleRangeCIs(likelyRangeSet,2)];

% Derivatives of all current data over injection range
injectionDataDiff = diff(currentVec(rangeCI(1):rangeCI(2), :), 1);

% All current data over injection range
injectionDataVals = currentVec(rangeCI(1):rangeCI(2), :);

% All current data outside injection range
noInjectionDataVals = currentVec(1:rangeCI(1), :);
noInjectionDataVals = [noInjectionDataVals; currentVec(rangeCI(2):end, :)];

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

isCI = false;
% Determine if CI:  1) Are the injection derivatives roughly the same?
%                   2) Are the injection derivatives roughly zero?
%                   3) Are the injections different from the rest of the recording?
%                   4) Is the injection range appropriately sized?
%                   2) Are there multiple sweeps being done?
if reduction < swpSpacingTolerance
    if injectionDataDiffMeans < zeroTolerance
        if diffInjection > zeroTolerance
            if rangeCI(2) - rangeCI(1) > minCIRangeTolerance
                if nSweeps > 1
                    isCI = true;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

OLD CODE:

if size(ciOptions,2) == 2                                                % if only 2 options remain, assume they represent start and end
       rangeCI = [ciOptions(1) ciOptions(2)];                        % may not be valid if max/max2/min/min2 did not find accurate points
elseif size(ciOptions,2) == 1 || size(ciOptions,2) == 0                % no valid CI range found
       rangeCI = [];
else
end

vcc = identify_channels(abfdata);
currentVec = squeeze(abfdata(:, find(vcc == 2), :));

% second derivative of current data for the max sweep (if CI, all CI should be on same time interval)
deriv2_cd = diff(currentVec(:,max_swp),2);

% [~, ci1] = max(deriv2_cd);

% [~, ci2] = max2(deriv2_cd);

% [~, ci3] = min(deriv2_cd);

% [~, ci4] = min2(deriv2_cd);
ciOptions = sort([ci1 ci2 ci3 ci4 ci5 ci6 ci7 ci8]);
% Store likely set as CI range
rangeCI = [possibleRangeCIs(likelyRangeSet,1) possibleRangeCIs(likelyRangeSet,2)];

% All current data over injection range
injection_data = currentVec(rangeCI(1):rangeCI(2), :);

% Get averages of each sweeps
avgs_byswp = mean(injection_data, 1);

% Difference between sweep averages should be minimal
reduction = abs(diff(avgs_byswp,2));

% Determine if CI:  1) Are the sweep averages roughly the same?
%                   2) Are there multiple sweeps being done?
if reduction < swpSpacingTolerance & size(abfdata,3) > 1
    isCI = true;
else
    isCI = false;
end

function [isCI, rangeCI] = identify_CI_protocol (abfdata, siUs)

%       /home/Matlab/Brians_Functions/identify_channels.m

% Get the channel types for each channel
channelTypes = identify_channels(abfdata);

% Find the index of the current in the .abf data
idxCurrent = find(strcmpi('Current', channelTypes));

% Extract the current vector
currentVec = squeeze(abfdata(:, idxCurrent, :));

% Apply lowpass filter
currentVec = freqfilter(currentVec, 50, siUs * 10^-6, 'FilterType', 'low');

if size(abfdata,3) > 1

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
if ~isdeployed
    addpath(fullfile(functionsdirectory, '/Adams_Functions/'));
                                            % for freqfilter.m
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%