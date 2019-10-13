function hinf = m3ha_compute_hinf_IT (v, shifth, slopeh)
% Compute the steady state value of the inactivation gating variable of IT
% Usage: hinf = m3ha_compute_hinf_IT (v, shifth, slopeh)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       cd/m3ha_compute_and_plot_IT.m
%
% File History:
% 2017-08-05 Created

% Parameters from Destexhe et al., 1998
vHalf = -81;                    % V_1/2 [mV] for minf
k = 4.0;                        % slope viewed sideways [mV] for minf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update vHalf & k
vHalf = vHalf + shifth;
k = k * slopeh;

% Compute minf
hinf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

