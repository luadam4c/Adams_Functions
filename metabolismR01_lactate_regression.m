% metabolismR01_lactate_regression.m
%% For interpolating lactate concentrations

% File History:
%   2020-02-25 Created by Adam Lu

%% Hard-coded parameters
directory = '/media/shareX/2020marchR01/lactate_regression';
standardCurveFileName = 'lactate_standard_curve.csv';
measurementsFileName = 'lactate_measurements.csv';
vStr = 'V';
lactateUmStr = 'lactateUm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find the full path to tables
[standardCurveFilePath, measurementsFilePath] = ...
    argfun(@(a) fullfile(directory, a), ...
            standardCurveFileName, measurementsFileName);

% Import tables
[standardCurveTable, measurementsTable] = ...
    argfun(@readtable, standardCurveFilePath, measurementsFilePath);

% Construct new sheetPath
sheetPath = replace(measurementsFilePath, '.csv', '_complete.csv');

%% Extract values
% Extract values from standard curve
voltStandards = standardCurveTable.(vStr);
lactateStandards = standardCurveTable.(lactateUmStr);

% Construct linear model
linearModel = fitlm(voltStandards, lactateStandards);

% Extract measurement values
voltValues = measurementsTable.(vStr);

% Compute corresponding voltages
lactateValues = predict(linearModel, voltValues);

% Add to table
measurementsTable = addvars(measurementsTable, lactateValues, 'NewVariableNames', lactateUmStr);

% Write table
writetable(measurementsTable, sheetPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%