% clc2_plot_measures.m
%% Plots all measures of interest across slices

% Requires:
%       cd/combine_variables_across_tables.m

% File History:
% 2019-03-15 Created by Adam Lu
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hard-coded parameters
% File patterns
sliceFilePattern = '.*slice.*';

% Must be consistent with parse_multiunit.m
varsToPlot = {'oscIndex1', 'oscIndex2', 'oscIndex3', 'oscIndex4', ...
                'oscPeriod1Ms', 'oscPeriod2Ms', ...
                'oscDurationSec', 'spikeCountTotal'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find all files with the pattern *slice*_params in the file name
[~, sliceParamSheets] = all_files('Keyword', 'slice', 'Suffix', 'params');

% Read all slice parameter spreadsheets
sliceParamsTables = cellfun(@readtable, sliceParamSheets, ...
                            'UniformOutput', false);

% Combine variables across tables
% measureTables = combine_variables_across_tables(sliceParamsTables, ...
%                                                 'VariableNames', varsToPlot);

%% Do the job
% TODO: Use plot_table in a different form: plot all columns together


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%