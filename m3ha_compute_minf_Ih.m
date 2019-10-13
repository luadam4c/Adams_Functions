function [minf, vHalf, k] = m3ha_compute_minf_Ih (v, shiftm)
% Compute the steady state value of the activation gating variable of Ih
% Usage: [minf, vHalf, k] = m3ha_compute_minf_Ih (v, shiftm)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       cd/m3ha_compute_and_plot_Ih.m
%
% File History:
% 2017-08-06 Created

% Parameters from Santoro et al., 2000
vHalf = -82;                    % V_1/2 [mV] for minf
k = 5.5;                        % slope viewed sideways [mV] for minf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update vHalf
vHalf = vHalf + shiftm;

% Compute minf
minf = boltzmann(v, vHalf, k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

