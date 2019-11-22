function current = compute_elcurr (conductance, voltage, reversalPotential)
%% Computes electrode current from conductance & voltage
% Usage: current = compute_elcurr (conductance, voltage, reversalPotential)
%
% Used by:    
%       cd/m3ha_trace_comparison.m

% File History:
% 2016-11-07 Moved from /media/adamX/m3ha/data_dclamp/trace_comparison.m
% 2018-07-09 Renamed variables to be more descriptive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computes the electrode current (positive inward)
%   Note: I = -g(V-E)
current = - conductance .* (voltage - reversalPotential);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%