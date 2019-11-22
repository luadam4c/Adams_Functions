function [vec1_rescaled, scaling_factor] = rescale_vec (vec1, vec2, precision)
% Rescale a vector (vec1) to be in the same ballpark as another vector (vec2),
% 	using the absolute maximum value for comparison
% Arguments:
%	vec1	  - vector to be rescaled
%	vec2	  - vector to compare to
%	precision - precision of scaling factor
%
% Used by:	
%		/media/adamX/m3ha/data_dclamp/trace_comparison.m
%		/media/adamX/m3ha/data_dclamp/m3ha_resave_sweeps.m
%
% 2016-11-07 Created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

absmax1 = max(abs(vec1));
absmax2 = max(abs(vec2));
scaling_factor = round((absmax2/absmax1) / precision) * precision;
vec1_rescaled = scaling_factor * vec1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
