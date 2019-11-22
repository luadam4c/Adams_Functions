function G = compute_gabab_conductance (t, Tlast, amp, Trise, TfallFast, TfallSlow, w)
%% Computes theoretical conductance curve for the GABA_B IPSC used by dynamic clamp
% Arguments:
% 		t	- Time vector (ms)
% 		Tlast	- Time of IPSC start (ms)
% 		amp	- Amplitude (nS)
% 		Trise	- Time constant for rising phase (ms)
% 		TfallFast - Time constant for fast falling phase (ms)
% 		TfallSlow - Time constant for slow falling phase (ms)
% 		w 	- weight of fast falling phase in the falling phase
%
% Used by:	
%		cd/m3ha_trace_comparison.m
%		/media/adamX/m3ha/data_dclamp/m3ha_resave_sweeps.m

% File History:
% 2016-11-07 Moved from /media/adamX/m3ha/data_dclamp/trace_comparison.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndp = length(t);		% number of time points
sims = t(2) - t(1);		% sampling interval (ms)
Tlast_ind = round(Tlast/sims);
G = zeros(ndp, 1);
for k = 1:ndp
	if k <= Tlast_ind
		G(k) = 0;
	else
		Ron(k) = exp(-(t(k)-Tlast)/Trise);
		RoffFast(k) = exp(-(t(k)-Tlast)/TfallFast);
		RoffSlow(k) = exp(-(t(k)-Tlast)/TfallSlow);
		G(k) = amp * (1 - Ron(k))^ 8 * (RoffFast(k)*w + RoffSlow(k)*(1-w));  
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
