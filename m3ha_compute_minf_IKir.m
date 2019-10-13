function minf = m3ha_compute_minf_IKir (v)
% Compute the steady state value of the activation gating variable of IKir
% Usage: minf = m3ha_compute_minf_IKir (v)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_IKir.m
%
% File History:
% 2017-08-06 Created

% Parameters from Amarillo et al., J Neurophysiol, 2014
vHalf = -97.9;                  % V_1/2 [mV] for minf
k = 9.7;                        % slope viewed sideways [mV] for minf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute minf
minf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
