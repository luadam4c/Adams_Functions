function m2inf = m3ha_compute_m2inf_IA (v)
% Compute the steady state value of the activation gating variable of IA
% Usage: m2inf = m3ha_compute_m2inf_IA (v)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       cd/m3ha_compute_and_plot_IA.m
%
% File History:
% 2017-08-06 Created

% Parameters from Huguenard & Prince, 1992
vHalf = -36;                    % V_1/2 [mV] for m2inf
k = -20;                        % slope viewed sideways [mV] for m2inf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute m2inf
m2inf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
