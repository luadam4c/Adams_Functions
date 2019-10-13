function tauh = m3ha_compute_tauh_IT (v, shifth, slopeh, celsius)
% Compute the time constant for the inactivation gating variable of IT
% Usage: tauh = m3ha_compute_tauh_IT (v, shifth, slopeh, celsius)
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compare_and_plot_IT.m
%
% File History:
% 2017-08-05 Created

% Parameters from Coulter et al., 1989
qh = 2.8;               % Q10 for inactivation

% Parameters from Destexhe et al., 1998
Trecord = 23;           % temperature of the voltage clamp experiments
theta = -80;            % threshold [mV]
v1 = -467;              % [mV]
k1 = 66.6;              % [mV]
tau2 = 28;              % [ms]
v2 = -22;               % [mV]
k2 = -10.5;             % [mV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute phih
phih = qh ^ ((celsius - Trecord)/10);

% Update theta, v1, v2, phih
theta = theta + shifth;
v1 = v1 + shifth; %% TODO: Does this make sense?
v2 = v2 + shifth;
% k1 = k1 * slopeh; %% TODO: Why not this?
% k2 = k2 * slopeh; %% TODO: Why not this?
phih = phih * slopeh;

% Compute tauh
tauh = zeros(size(v));
for i = 1:length(v)
    if  v(i) < theta
        tauh(i) = exp((v(i)-v1)./k1)./phih;
    else
        tauh(i) = (tau2 + exp((v(i)-v2)./k2))./phih;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

