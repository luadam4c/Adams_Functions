function tauh2 = m3ha_compute_tauh2_IA (v, celsius)
% Compute the time constant for the inactivation gating variable of IA
% Usage: tauh2 = m3ha_compute_tauh2_IA (v, celsius)
%
% Used by:    
%       cd/m3ha_compute_and_plot_IA.m
%
% File History:
% 2017-08-06 Created

% Parameters from Huguenard et al., 1991
qh = 2.8;               % Q10 for inactivation

% Parameters from Huguenard & Prince, 1992
Trecord = 23;           % temperature of the voltage clamp experiments
theta = -73;            % threshold [mV]
v1 = -46;               % [mV]
k1 = 5.0;               % [mV]
v2 = -238;              % [mV]
k2 = -37.5;             % [mV]
tau3 = 60;              % [ms]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute phih
phih = qh ^ ((celsius - Trecord)/10);

% Compute tauh2
tauh2 = zeros(size(v));
for i = 1:length(v)
    if  v(i) < theta
        tauh2(i) = 1/(exp((v(i)-v1)/k1) + exp((v(i)-v2)/k2))/phih;
    else
        tauh2(i) = tau3/phih;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

