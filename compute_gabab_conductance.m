function gGabab = compute_gabab_conductance (tVec, timeStart, amplitude, ...
                                    tauRise, tauFallFast, tauFallSlow, weight)
%% Computes a the conductance over time for a GABAB-IPSC
% Usage: gGabab = compute_gabab_conductance (tVec, timeStart, amplitude, ...
%                                   tauRise, tauFallFast, tauFallSlow, weight)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       gGabab      - GABAB-IPSC conductance (nS)
%                   specified as a numeric vector
%
% Arguments:
%       tVec        - time vector (ms)
%                   must be a numeric vector
%       timeStart   - time of IPSC start (ms)
%                   must be a numeric scalar
%       amplitude   - amplitude (nS)
%                   must be a numeric scalar
%       tauRise     - time constant for rising phase (ms)
%                   must be a numeric scalar
%       tauFallFast - time constant for fast falling phase (ms)
%                   must be a numeric scalar
%       tauFallSlow - time constant for slow falling phase (ms)
%                   must be a numeric scalar
%       weight      - weight of fast falling phase in the falling phase
%                   must be a numeric scalar
%
% Requires:
%       cd/find_closest.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%
% Used by:    
%       cd/m3ha_resave_sweeps.m
%       cd/m3ha_trace_comparison.m

% File History:
% 2016-11-07 Moved from /media/adamX/m3ha/data_dclamp/trace_comparison.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Make sure time vectors are column vectors
tVec = force_column_vector(tVec);

% Make sure parameters are row vectors
[timeStart, amplitude, tauRise, tauFallFast, tauFallSlow, weight] = ...
    argfun(@force_row_vector, ...
            timeStart, amplitude, tauRise, tauFallFast, tauFallSlow, weight);

% Compute the shifted time vectors
timeShifted = tVec - timeStart;

% Compute the exponential components
[expRise, expFallFast, expFallSlow] = ...
    argfun(@(x) exp(-timeShifted ./ x), tauRise, tauFallFast, tauFallSlow);

% Compute the GABA-B conductance
gGabab = amplitude .* (1 - expRise) .^ 8 * ...
                (expFallFast .* weight + expFallSlow .* (1 - weight));  

% Compute the index of the starting time
idxTimeStart = find_closest(tVec, timeStart);

% For each column, set everything before timeStart to be zero
for iColumn = 1:size(gGabab, 2)
    gGabab(1:idxTimeStart(iColumn), iColumn) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%