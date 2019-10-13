function minf = m3ha_compute_minf_IT (v, shiftm, slopem)
% Compute the steady state value of the activation gating variable of IT
% Usage: minf = m3ha_compute_minf_IT (v, shiftm, slopem)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       cd/m3ha_compute_and_plot_IT.m
%
% File History:
% 2017-08-05 Created

% Parameters from Destexhe et al., 1998
vHalf = -57;                    % V_1/2 [mV] for minf
k = -6.2;                       % slope viewed sideways [mV] for minf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update vHalf & k
vHalf = vHalf + shiftm;
k = k * slopem;

% Compute minf
minf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

