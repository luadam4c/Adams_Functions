function hinf = m3ha_compute_hinf_IA (v)
% Compute the steady state value of the inactivation gating variable of IA
% Usage: hinf = m3ha_compute_hinf_IA (v)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_IA.m
%
% File History:
% 2017-08-06 Created

% Parameters from Huguenard & Prince, 1992
vHalf = -78;                    % V_1/2 [mV] for hinf
k = 6;                          % slope viewed sideways [mV] for hinf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute hinf
hinf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
