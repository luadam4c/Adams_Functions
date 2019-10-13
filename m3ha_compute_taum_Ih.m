function taum = m3ha_compute_taum_Ih (v, shiftm, celsius)
% Compute the time constant for the activation gating variable of IT
% Usage: taum = m3ha_compute_taum_Ih (v, shiftm, celsius)
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_IT.m
%
% File History:
% 2017-08-05 Created

% Parameters from Santoro & Tibbs, 1999
qm = 4.0;               % Q10 for activation

% Parameters from Destexhe et al., 1998
Trecord = 34;           % temperature of the voltage clamp experiments
tau0 = 8e-4;            % [ms]
tau1 = 3.5e-6;          % [ms]
v1 = 0;                 % [mV]
r1 = -0.05787;          % [/mV]
lambda2 = -1.87;        % [1]
r2 = 0.0701;            % [/mV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute phim
phim = qm ^ ((celsius - Trecord)/10);

% Update v1
v1 = v1 + shiftm;       %% TODO: Does this make sense?

% Compute taum
taum = (1./(tau0 + tau1 * exp(r1 * (v - v1)) ...
        + exp(lambda2 + r2 * (v - v1)))) ./ phim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

