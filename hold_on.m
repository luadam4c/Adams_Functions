function wasHold = hold_on
%% Holds on and returns previous status
% Usage: wasHold = hold_on
% Explanation:
%       TODO
%
% Example(s):
%       wasHold = hold_on;
%       hold_off(wasHold);
%
% Outputs:
%       wasHold     - whether the current axes was held on
%                   specified as a logical scalar
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m
%       cd/plot_autocorrelogram.m
%       cd/plot_ball_stick.m
%       cd/plot_chevron.m
%       cd/plot_grouped_jitter.m
%       cd/plot_raw_multiunit.m
%       cd/plot_spike_density_multiunit.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-10-02 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
if ~ishold
    wasHold = false;
    hold on
else
    wasHold = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%