function hinf = m3ha_compute_hinf_INaP (v)
% Compute the steady state value of the inactivation gating variable of INaP
% Usage: hinf = m3ha_compute_hinf_INaP (v)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_INaP.m
%
% File History:
% 2017-08-06 Created

% Parameters from Wu et al., 2005
vHalf = -58.7;                  % V_1/2 [mV] for hinf
k = 14.2;                       % slope viewed sideways [mV] for hinf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute hinf
hinf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
