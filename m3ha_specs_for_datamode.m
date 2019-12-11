function [fileSuffix, titleModifier] = m3ha_specs_for_datamode (dataMode)
%% Specifications depending on dataMode for the GAT blockade project
% Usage: [fileSuffix, titleModifier] = m3ha_specs_for_datamode (dataMode)
%
% Used by:
%       cd/m3ha_compute_and_plot_statistics.m
%       cd/m3ha_estimate_passive_params.m
%       cd/m3ha_plot_histograms_refine_threshold.m
%       cd/m3ha_plot_correlations.m

% File History:
% 2016-10-31 Created
% 2017-05-22 AL - Renamed function SpecsForFitmode() -> specs_for_fitmode()
% 2018-10-04 AL - Moved to Adams_Functions
% 2018-10-04 AL - Renamed function specs_for_fitmode -> m3ha_specs_for_datamode
% 2018-10-04 AL - Added dataMode == -1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dataMode == -1
    fileSuffix = '';
    titleModifier = '';
elseif dataMode == 0
    fileSuffix = '_all';
    titleModifier = '(all)';
elseif dataMode == 1
    fileSuffix = '_100-400all';
    titleModifier = '(100%, 200%, 400% g incr)';
elseif dataMode == 2
    fileSuffix = '_tofit';
    titleModifier = '(for fitting)';
else
    error('The data mode %d is unrecognized!!\n', dataMode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       cd/find_passive_params.m

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
