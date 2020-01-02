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
%
% Used by:    
%       cd/m3ha_resave_sweeps.m
%       cd/m3ha_trace_comparison.m

% File History:
% 2016-11-07 Moved from /media/adamX/m3ha/data_dclamp/trace_comparison.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Compute the number of samples
nSamples = numel(tVec);

% Compute the index of the starting time
idxTimeStart = find_closest(tVec, timeStart);

gGabab = zeros(nSamples, 1);
for k = 1:nSamples
    if k <= idxTimeStart
        gGabab(k) = 0;
    else
        Ron(k) = exp(-(tVec(k)-timeStart)/tauRise);
        RoffFast(k) = exp(-(tVec(k)-timeStart)/tauFallFast);
        RoffSlow(k) = exp(-(tVec(k)-timeStart)/tauFallSlow);
        gGabab(k) = amplitude * (1 - Ron(k))^ 8 * (RoffFast(k)*weight + RoffSlow(k)*(1-weight));  
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

siMs = tVec(2) - tVec(1);        % sampling interval (ms)
idxTimeStart = round(timeStart/siMs);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%