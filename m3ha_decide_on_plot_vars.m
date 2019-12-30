function [varsToPlot, varLabelsToPlot, yLimitsToPlot, varIsLogToPlot] = ...
                m3ha_decide_on_plot_vars
%% Decides on the error and parameters to plot
% Usage: [varsToPlot, varLabelsToPlot, yLimitsToPlot, varIsLogToPlot] = ...
%               m3ha_decide_on_plot_vars
% Explanation:
%       TODO
%
% Example(s):
%       [vars, labels, yLims, isLog] = m3ha_decide_on_plot_vars
%
% Outputs:
%       TODO
%
% Arguments:
%
% Requires:
%       cd/find_first_match.m
%       cd/force_column_cell.m
%       cd/m3ha_neuron_create_default_params.m
%
% Used by:
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_rank_neurons.m

% File History:
% 2019-12-30 Moved from m3ha_neuron_choose_best_params.m
% 

%% Hard-coded parameters
errorsToPlot = {'totalError'; 'avgSwpError'; 'ltsMatchError'; ...
                'avgLtsAmpError'; 'avgLtsDelayError'; 'avgLtsSlopeError'};
errorLabelsToPlot = {'Total Error'; 'Sweep Error'; 'Match Error'; ...
            'Amp Error'; 'Time Error'; 'Slope Error'};
errorYLimits = [0; Inf];
paramsToPlot = { 'diamSoma'; 'LDend'; 'diamDend'; 'gpas'; 'epas'; ...
                'pcabarITSoma'; 'ghbarIhSoma'; 'gkbarIKirSoma'; ...
                'gkbarIASoma'; 'gnabarINaPSoma'; ...
                'pcabarITDend1'; 'ghbarIhDend1'; 'gkbarIKirDend1'; ...
                'gkbarIADend1'; 'gnabarINaPDend1'; ...
                'pcabarITDend2'; 'ghbarIhDend2'; 'gkbarIKirDend2'; ...
                'gkbarIADend2'; 'gnabarINaPDend2'};
paramLabelsToPlot = paramsToPlot;

% Note: must be consistent with m3ha_neuron_create_default_params.m
lowerBoundStr = 'LowerBound';
upperBoundStr = 'UpperBound';
isLogStr = 'IsLog';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with NEURON parameters
% Create the default NEURON parameters table
defaultTable = m3ha_neuron_create_default_params;

% Extract information from the table
neuronParamNames = defaultTable.Properties.RowNames;
neuronParamsLowerBound = defaultTable.(lowerBoundStr);
neuronParamsUpperBound = defaultTable.(upperBoundStr);
neuronParamsIsLog = defaultTable.(isLogStr);

% Decide on NEURON parameter y axis limits
neuronParamsYLimits = transpose([neuronParamsLowerBound, neuronParamsUpperBound]);

% Get the original index for each parameter to plot
indParamsToPlot = find_first_match(paramsToPlot, neuronParamNames, ...
                    'MatchMode', 'exact', 'IgnoreCase', true);

% Extract the y limits for each parameter
paramYLimits = force_column_cell(neuronParamsYLimits(:, indParamsToPlot));

% Determine whether each parameter should be plotted in a log scale
paramIsLog = neuronParamsIsLog(indParamsToPlot);

%% Combine errors and parameters
% Decide on the errors and parameters to plot
varsToPlot = vertcat(errorsToPlot(2:6), paramsToPlot);
varLabelsToPlot = vertcat(errorLabelsToPlot(2:6), paramLabelsToPlot);

% Construct all error and parameter y limits
yLimitsToPlot = [repmat({errorYLimits}, 5, 1); paramYLimits];
varIsLogToPlot = [false(5, 1); paramIsLog];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
