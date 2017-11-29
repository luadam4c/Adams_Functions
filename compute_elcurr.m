function Curr = compute_elcurr (G, V, E_rev)
% Compute electrode current from conductance & voltage
%
% Used by:	
%		/media/adamX/m3ha/data_dclamp/trace_comparison.m
%
% 2016-11-07 Moved from /media/adamX/m3ha/data_dclamp/trace_comparison.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Curr = - G .* (V - E_rev);

