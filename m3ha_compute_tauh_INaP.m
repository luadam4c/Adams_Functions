function tauh = m3ha_compute_tauh_INaP (v, celsius)
% Compute the time constant for the inactivation gating variable of INaP
% Usage: tauh = m3ha_compute_tauh_INaP (v, celsius)
%
% Requires:
%       /home/Matlab/boltzmann.m
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_INaP.m
%
% File History:
% 2017-08-06 Created

% Parameters assumed by Amarillo et al., 2014
qh = 3;               % Q10 for inactivation

% Parameters from Wu et al., 2005
Trecord = 23;           % temperature of the voltage clamp experiments
tau0 = 1000;            % [ms]
tau1 = 10000;           % [ms]
vHalf = -60;            % [mV]
k = 10;                 % [mV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute phih
phih = qh ^ ((celsius - Trecord)/10);

% Compute tauh
tauh = (tau0 + (tau1 .* boltzmann(v, vHalf, k) ) ) ./ phih;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

