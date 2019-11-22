function [fileSuffix, titleModifier] = m3ha_specs_for_fitmode (fitMode)
%% Specifications depending on fitMode for the GAT blockade project
% Usage: [fileSuffix, titleModifier] = m3ha_specs_for_fitmode (fitMode)
%
% Used by:
%       cd/m3ha_dclampdatalog_analyze.m
%       cd/m3ha_dclampPassiveFitter.m
%       cd/m3ha_PlotHistogramsRefineThreshold.m
%       cd/m3ha_PlotCorrelations.m

% File History:
% 2016-10-31 Created
% 2017-05-22 AL - Renamed function SpecsForFitmode() -> specs_for_fitmode()
% 2018-10-04 AL - Moved to Adams_Functions
% 2018-10-04 AL - Renamed function specs_for_fitmode -> m3ha_specs_for_fitmode
% 2018-10-04 AL - Added fitMode == -1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fitMode == -1
    fileSuffix = '';
    titleModifier = '';
elseif fitMode == 0
    fileSuffix = '_all';
    titleModifier = '(all)';
elseif fitMode == 1
    fileSuffix = '_100-400all';
    titleModifier = '(100%, 200%, 400% g incr)';
elseif fitMode == 2
    fileSuffix = '_tofit';
    titleModifier = '(for fitting)';
else
    error('The fit mode %d is unrecognized!!\n', fitMode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       cd/find_passive_params.m

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
