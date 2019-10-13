function m1inf = m3ha_compute_m1inf_IA (v)
% Compute the steady state value of the activation gating variable of IA
% Usage: m1inf = m3ha_compute_m1inf_IA (v)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       cd/m3ha_compute_and_plot_IA.m
%
% File History:
% 2017-08-06 Created

% Parameters from Huguenard & Prince, 1992
vHalf = -60;                    % V_1/2 [mV] for m1inf
k = -8.5;                       % slope viewed sideways [mV] for m1inf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute m1inf
m1inf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
