function taum = m3ha_compute_taum_IA (v, celsius)
% Compute the time constant for the activation gating variable of IA
% Usage: taum = m3ha_compute_taum_IA (v, celsius)
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_IA.m
%
% File History:
% 2017-08-06 Created

% Parameters from Huguenard et al., 1991
qm = 2.8;               % Q10 for activation

% Parameters from Huguenard & Prince, 1992
Trecord = 23;           % temperature of the voltage clamp experiments
tau0 = 0.37;            % [ms]
v1 = -35.8;             % [mV]
v2 = -79.7;             % [mV]
k1 = 19.7;              % [mV]
k2 = -12.7 ;            % [mV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute phim
phim = qm ^ ((celsius - Trecord)/10);

% Compute taum
taum = (tau0 + 1./(exp((v-v1)./k1)+exp((v-v2)./k2))) ./ phim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

