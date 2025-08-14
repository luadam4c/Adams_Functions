function tpToPrevTrue = compute_timepoints_to_previous(isTrue)
%% Returns a vector with the number of time points to previous TRUE index, for each TRUE index
% Usage: tpToPrevTrue = compute_timepoints_to_previous(isTrue)
% Explanation:
%       Computes the number of time points to previous TRUE index,
%           for each TRUE index in a logical vector
%
% Example(s):
%       tpToPrevTrue = compute_timepoints_to_previous([0, 1, 1, 0, 1, 1])
%
% Outputs:
%       tpToPrevTrue    - number of time points to previous TRUE index, 
%                           for each TRUE index
%                       specified as a numeric vector
%
% Arguments:    
%       isTrue      - a vector of whether each index (time point) is TRUE
%                   must be a logical vector
%       
% Requires:
%
% Used by:
%       \Shared\scAAV\analyze_reachr_motion.m

% File History:
% 2025-08-13 Modified from getPrevStim.m by Jeff Moore

% Find the TRUE indices
indStim = find(isTrue);

% Create a vector of NaNs with the same size as the TRUE indices
tpToPrevTrue = nan(size(indStim));

% The first index has no previous index
tpToPrevTrue(1) = indStim(1);

% Count the number of indices to the previous index
tpToPrevTrue(2:end) = indStim(2:end) - indStim(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%