function minf = m3ha_compute_minf_INaP (v)
% Compute the steady state value of the activation gating variable of INaP
% Usage: minf = m3ha_compute_minf_INaP (v)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       cd/m3ha_compute_and_plot_INaP.m
%
% File History:
% 2017-08-06 Created

% Parameters from Wu et al., 2005
vHalf = -57.9;                  % V_1/2 [mV] for minf
k = -6.4;                       % slope viewed sideways [mV] for minf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute minf
minf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
