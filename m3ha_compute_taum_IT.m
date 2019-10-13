function taum = m3ha_compute_taum_IT (v, shiftm, slopem, celsius)
% Compute the time constant for the activation gating variable of IT
% Usage: taum = m3ha_compute_taum_IT (v, shiftm, slopem, celsius)
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_IT.m
%
% File History:
% 2017-08-05 Created

% Parameters from Coulter et al., 1989
qm = 3.6;               % Q10 for activation

% Parameters from Destexhe et al., 1998
Trecord = 23;           % temperature of the voltage clamp experiments
tau0 = 0.612;           % [ms]
v1 = -132;              % [mV]
v2 = -16.8;             % [mV]
k1 = -16.7;             % [mV]
k2 = 18.2;              % [mV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute phim
phim = qm ^ ((celsius - Trecord)/10);

% Update v1, v2, phim
v1 = v1 + shiftm; %% TODO: Does this make sense?
v2 = v2 + shiftm;
% k1 = k1 * slopem; %% TODO: Why not this?
% k2 = k2 * slopem; %% TODO: Why not this?
phim = phim * slopem;

% Compute taum
taum = (tau0 + 1./(exp((v-v1)./k1)+exp((v-v2)./k2))) ./ phim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

