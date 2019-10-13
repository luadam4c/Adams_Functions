function p = boltzmann(v, vHalf, k)
%% Computes the sigmoidal Boltzmann function
% Usage: p = boltzmann(v, vHalf, k)
%
% Outputs:
%       p       - the ordinate values
% Arguments:
%       v       - the abscissa values
%               must be a numeric array
%       vHalf   - half-maximum point
%               must be a numeric scalar
%       k       - slope when viewed sideways
%               must be a numeric scalar
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_minf_IT.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_hinf_IT.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_minf_Ih.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_minf_IKir.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_m1inf_IA.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_m2inf_IA.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_hinf_IA.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_minf_INaP.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_hinf_INaP.m
%       /media/adamX/m3ha/optimizer4gabab/m3ha_compute_tauh_INaP.m
% 
% File History:
% 2017-08-06 Created
% TODO: Input parser
% 

p = 1 ./ ( 1 + exp( (v - vHalf) ./ k ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
