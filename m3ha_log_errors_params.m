function m3ha_log_errors_params (logFileName, outparams, err, simplexOut)
%% Log errors and parameter values
% Usage: m3ha_log_errors_params (logFileName, outparams, err, simplexOut)
%
% Requires:
%       cd/save_params.m
%
% Used by:
%       cd/m3ha_fminsearch3.m
%       ~/m3ha/optimizer4gabab/m3ha_optimizer_4compgabab.m

% File History:
% 2017-01-17 Created
% 2017-01-25 Added concise version
% 2017-04-21 Added outparams.fitreg and outparams.fitregCpr
% 2017-05-15 Added currentparams to save to .p file
% 2017-05-16 Replaced tolx and tolf with relativeParamTolerance and 
%               relativeErrorTolerance
% 2017-05-17 Modified what to show and added outparams.logSwpErrFlag, 
%        outparams.logSwpWeightsFlag, outparams.logLtsWeightsFlag, 
%        outparams.logSwpEndsFlag
% 2017-05-17 If lts not computed, do not show errors
% 2018-08-10 Updated fitreg -> fitWindow
% 2018-11-16 Reorganized code and made logFileName the first argument
% TODO: Make each iteration a new column instead of a new row

% Fields for lts errors to check if exist in err
ltsErrorFields = {'avgLtsAmpError', 'avgLtsDelayError', 'avgLtsSlopeError'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Check if logging simplex version
if nargin < 4 || isempty(simplexOut)
    isSimplex = false;
else
    isSimplex = true;
end

% Decide on the neuron parameters to extract
if isSimplex
    neuronParamsTable = simplexOut.neuronParamsTable;
else 
    neuronParamsTable = outparams.neuronParamsTable;
end

% Extract info from outparams
simMode = outparams.simMode;
logSwpErrFlag = outparams.logSwpErrFlag;
logSwpWeightsFlag = outparams.logSwpWeightsFlag;
logLtsWeightsFlag = outparams.logLtsWeightsFlag;
logSwpEndsFlag = outparams.logSwpEndsFlag;
saveParamsFlag = outparams.saveParamsFlag;
outFolderName = outparams.outFolder;
prefix = outparams.prefix;
simplexNum = outparams.simplexNum;
fitWindowCpr = outparams.fitWindowCpr;
fitWindowIpscr = outparams.fitWindowIpscr;
sweepWeightsIpscr = outparams.sweepWeightsIpscr;
ltsFeatureWeights = outparams.ltsFeatureWeights;
lts2SweepErrorRatio = outparams.lts2SweepErrorRatio;

% Extract info from err
swpErrors = err.swpErrors;
avgSwpError = err.avgSwpError;
ltsAmpErrors = err.ltsAmpErrors;
ltsDelayErrors = err.ltsDelayErrors;
ltsSlopeErrors = err.ltsSlopeErrors;
avgLtsAmpError = err.avgLtsAmpError;
avgLtsDelayError = err.avgLtsDelayError;
avgLtsSlopeError = err.avgLtsSlopeError;
avgLtsError = err.avgLtsError;
totalError = err.totalError;

% Extract other information from simplexOut
if isSimplex
    ctIterations = simplexOut.ctIterations;
    maxIterations = simplexOut.maxIterations;
    lastHow = simplexOut.lastHow;
    ctEvals = simplexOut.ctEvals;
    maxFunctionEvaluations = simplexOut.maxFunctionEvaluations;
    totalError = simplexOut.totalError;
    maxErrorChange = simplexOut.maxErrorChange;
    relativeErrorTolerance = simplexOut.relativeErrorTolerance;
    maxParamChange = simplexOut.maxParamChange;
    relativeParamTolerance = simplexOut.relativeParamTolerance;

    % Get the names and values of parameters in use
    paramsInUseValues = simplexOut.paramsInUseValues;
    paramsInUseNames = simplexOut.paramsInUseNames;

    % Count the number of parameters in use
    nInUse = length(paramsInUseValues);    
end

% Extract from neuronParamsTable
neuronParamNames = neuronParamsTable.Properties.RowNames;
neuronParamValues = neuronParamsTable.Value;
neuronParamLowerBounds = neuronParamsTable.LowerBound;
neuronParamUpperBounds = neuronParamsTable.UpperBound;

% Count the number of parameters
nParams = height(neuronParamsTable);

% Count the number of sweeps
nSweeps = length(sweepWeightsIpscr);

% Check if LTS error exists
hasLtsError = isfield(err, ltsErrorFields);

% Create a different file name for the concise version
conciseFileName = strrep(logFileName, '.csv', '_concise.csv');

% Decide on the output parameters spreadsheet name
if isSimplex
    sheetName = fullfile(outFolderName, [prefix, ...
                        '_simplexrun_', num2str(simplexNum), ...
                        '_iter_', num2str(ctIterations), '_params.csv']);
else 
    sheetName = fullfile(outFolderName, [prefix, '_params.csv']);
end


%% Do the job
% Create file if not exist and print log header
logFilePath = fullfile(outFolderName, logFileName);
if ~isfile(logFilePath)
    fid = fopen(logFilePath, 'w');        % create log file for writing
    fprintf(fid, '%s, ', 'Simulation Mode');
    if isSimplex
        fprintf(fid, '%s, %s, %s, %s, %s, %s, ', ...
            'ctIterations', 'maxIterations', 'how', ...
            'ctEvals', 'maxFunctionEvaluations', 'fval');
        fprintf(fid, '%s, %s, %s, %s, ', ...
            'Maximum relative error change', ...
            'Error tolerance (relative)', ...
            'Maximum relative parameter change', ...
            'Parameter change tolerance (relative)');
    end
    fprintf(fid, '%s, ', '# of sweeps');
    if isSimplex
        fprintf(fid, '%s, ', '# of parameters in use');
    end
    fprintf(fid, '%s, %s, ', ...
            'Total error', 'Normalized average sweep error');
    if all(hasLtsError)
        fprintf(fid, '%s, %s, %s, %s, ', ...
            'Normalized average LTS amp error', ...
            'Normalized average LTS time error', ...
            'Normalized average LTS slope error', ...
            'Normalized average LTS error');
    end
    fprintf(fid, repmat('%s, ', 1, nParams), neuronParamNames{:});
    if logSwpErrFlag
        strings = repmat({'sweep error (mV)'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
    end
    if all(hasLtsError) && logSwpErrFlag
        strings = repmat({'LTS amp error'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
        strings = repmat({'LTS time error'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
        strings = repmat({'LTS slope error'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
    end
    if logSwpWeightsFlag
        strings = repmat({'sweep weight'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
    end
    if all(hasLtsError) && logLtsWeightsFlag
        fprintf(fid, repmat('%s, ', 1, 3), ...
            'LTS amp weight', 'LTS time weight', ...
            'LTS slope weight');
        fprintf(fid, '%s, ', 'LTS to sweep error ratio');
    end
    fprintf(fid, '%s, %s, %s, %s, ', ...
        'cpr start', 'cpr end', 'sweep start', 'sweep end');
    %{
    if logSwpEndsFlag
        strings = repmat({'cpr start'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
        strings = repmat({'cpr end'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
        strings = repmat({'sweep start'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
        strings = repmat({'sweep end'}, 1, nSweeps);
        fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
    end
    %}
    fprintf(fid, '\n');
    fclose(fid);
end

% Log errors and params
fid = fopen(logFilePath, 'a');            % append
fprintf(fid, '%s, ', simMode);
if isSimplex
    fprintf(fid, '%6.4g, %6.4g, %s, %6.4g, %6.4g, %6.4g, ', ...
            ctIterations, maxIterations, lastHow, ctEvals, ...
            maxFunctionEvaluations, totalError);
    fprintf(fid, '%6.4g, %6.4g, %6.4g, %6.4g, ', ...
            maxErrorChange, relativeErrorTolerance, ...
            maxParamChange, relativeParamTolerance);
end
fprintf(fid, '%2.2g, ', nSweeps);
if isSimplex
    fprintf(fid, '%2.2g, ', nInUse);
end
fprintf(fid, '%6.4g, %6.4g, ', totalError, avgSwpError);
if all(hasLtsError)
    fprintf(fid, '%6.4g, %6.4g, %6.4g, %6.4g, ', ...
            avgLtsAmpError, avgLtsDelayError, avgLtsSlopeError, avgLtsError);
end
fprintf(fid, repmat('%6.4g, ', 1, nParams), neuronParamValues);
if logSwpErrFlag
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), swpErrors);
end
if all(hasLtsError) && logSwpErrFlag
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), ltsAmpErrors);
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), ltsDelayErrors);
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), ltsSlopeErrors);
end
if logSwpWeightsFlag
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), sweepWeightsIpscr);
end
if all(hasLtsError) && logLtsWeightsFlag
    fprintf(fid, repmat('%6.4g, ', 1, 3), ltsFeatureWeights);
    fprintf(fid, '%6.4g, ', lts2SweepErrorRatio);
end
fprintf(fid, '%6.4g, ', fitWindowCpr(1));
fprintf(fid, '%6.4g, ', fitWindowCpr(2));
fprintf(fid, '%6.4g, ', fitWindowIpscr(1));
fprintf(fid, '%6.4g, ', fitWindowIpscr(2));
%{
if logSwpEndsFlag
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), fitWindowCpr(:, 1)');
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), fitWindowCpr(:, 2)');
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), fitWindowIpscr(:, 1)');
    fprintf(fid, repmat('%6.4g, ', 1, nSweeps), fitWindowIpscr(:, 2)');
end
%}
fprintf(fid, '\n');
fclose(fid);

% Create concise version if not exist and print log header
conciseFilePath = fullfile(outFolderName, conciseFileName);
if ~isfile(conciseFilePath)
    fid = fopen(conciseFilePath, 'w');        % create
    if isSimplex
        fprintf(fid, '%s, %s, %s, %s, ', ...
            'ctIterations', 'how', 'ctEvals', 'Total error');
        fprintf(fid, '%s, %s, ', ...
            'Maximum error change', 'Maximum parameter change');
    else
        fprintf(fid, '%s, ', 'Total error');
    end
    fprintf(fid, '%s, ', 'Normalized average sweep error');
    if all(hasLtsError)
        fprintf(fid, '%s, %s, %s, %s, ', ...
            'Normalized average LTS amp error', ...
            'Normalized average LTS time error', ...
            'Normalized average LTS slope error', ...
            'Normalized average LTS error');
    end
    if isSimplex
        fprintf(fid, repmat('%s, ', 1, nInUse), paramsInUseNames{:});
    else
        fprintf(fid, repmat('%s, ', 1, nParams), neuronParamNames{:});
    end
    if all(hasLtsError) && logLtsWeightsFlag
        fprintf(fid, repmat('%s, ', 1, 3), ...
            'LTS amp weight', 'LTS time weight', ...
            'LTS slope weight');
        fprintf(fid, '%s, ', 'LTS to sweep error ratio');
    end
    fprintf(fid, '\n');
    fclose(fid);
end

% Log errors and params (concise version)
fid = fopen(conciseFilePath, 'a');        % append
if isSimplex
    fprintf(fid, '%6.4g, %s, %6.4g, %6.4g, ', ...
            ctIterations, lastHow, ctEvals, totalError);
    fprintf(fid, '%6.4g, %6.4g, ', ...
            maxErrorChange, maxParamChange);
else
    fprintf(fid, '%6.4g, ', totalError);
end
fprintf(fid, '%6.4g, ', avgSwpError);
if all(hasLtsError)
    fprintf(fid, '%6.4g, %6.4g, %6.4g, %6.4g, ', ...
            avgLtsAmpError, avgLtsDelayError, avgLtsSlopeError, avgLtsError);
end
if isSimplex
    fprintf(fid, repmat('%6.4g, ', 1, nInUse), paramsInUseValues);
else
    fprintf(fid, repmat('%6.4g, ', 1, nParams), neuronParamValues);
end
if all(hasLtsError) && logLtsWeightsFlag
    fprintf(fid, repmat('%6.4g, ', 1, 3), ltsFeatureWeights);
    fprintf(fid, '%6.4g, ', lts2SweepErrorRatio);
end
fprintf(fid, '\n');
fclose(fid);

%% Save parameters into a p file
%   Note: Always do this if not in simplex mode
if ~isSimplex || saveParamsFlag
    save_params(neuronParamsTable, 'FileName', sheetName);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE

fprintf(fid, '%2.2g, ', nSweeps);
fprintf(fid, '%6.4f, ', outparams.lts2SweepErrorRatio);
fprintf(fid, repmat('%6.4e, ', 1, numel(outparams.neuronParamValues)), outparams.neuronParamValues);
fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.swpw)), outparams.swpw);
fprintf(fid, '%6.4e, ', err.totalError);
fprintf(fid, repmat('%6.4e, ', 1, numel(err.swpErrors)), err.swpErrors);
fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.ltsFeatureWeights)), outparams.ltsFeatureWeights);
fprintf(fid, repmat('%6.4e, ', 1, 3), err.ltserrv, err.ltserrt, err.ltserrdvdtv);
fprintf(fid,'\n');

strings = repmat({'sweep left edge'}, 1, nSweeps);
fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
strings = repmat({'sweep right edge'}, 1, nSweeps);
fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});

fprintf(fid, repmat('%6.4g, ', 1, nSweeps), outparams.swpedges(:,1)');
fprintf(fid, repmat('%6.4g, ', 1, nSweeps), outparams.swpedges(:,2)');

fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.fitreg(:, 1))), outparams.fitreg(:, 1)');
fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.fitreg(:, 2))), outparams.fitreg(:, 2)');

if outparams.cprflag

if isSimplex
    fprintf(fid, repmat('%6.4g, ', 1, nParams), simplexOut.neuronParamValues);
else
    fprintf(fid, repmat('%6.4g, ', 1, nParams), outparams.neuronParamValues);
end

tocheck = {'ltserr_maxamp_val', 'ltserr_maxamp_time', 'ltserr_maxdiffamp_val'};
fprintf(fid, '%s, %s, %s', 'LTS amp error', 'LTS time error', 'LTS slope error');
fprintf(fid, '%6.4g, %6.4g, %6.4g, ', err.ltserr_maxamp_val, err.ltserr_maxamp_time, err.ltserr_maxdiffamp_val);

strings = repmat({'sweep error'}, 1, nSweeps);
fprintf(fid, repmat('%s, ', 1, nSweeps), strings{:});
fprintf(fid, repmat('%6.4g, ', 1, nSweeps), err.swpErrors);

for k = 1:numel(tocheck)
    if ~isfield(err, ltsErrorFields{k})
        err.(ltsErrorFields{k}) = NaN;
    end
end

hasLtsError = zeros(size(ltsErrorFields));
for k = 1:numel(ltsErrorFields)
    if isfield(err, ltsErrorFields{k})
        hasLtsError(k) = 1;
    end
end

pfilename = fullfile(outparams.outFolder, [outparams.prefix, '.p']);
%% Save parameters into a p file
currentparams = [neuronParamNames', num2cell(neuronParamValues'), ...
    num2cell(neuronParamLowerBounds'), num2cell(neuronParamUpperBounds')];
save(pfilename, 'currentparams');

pfilename = fullfile(outparams.outFolder, ...
    [outparams.prefix, ...
    '_simplexrun_', num2str(outparams.simplexNum), ...
    '_iter_', num2str(simplexOut.ctIterations), '.p']);

currentparams = [neuronParamNames', num2cell(neuronParamValues'), ...
    num2cell(neuronParamLowerBounds'), num2cell(neuronParamUpperBounds')];
save(pfilename, 'currentparams');

save_params(sheetName, neuronParamNames, neuronParamValues, ...
            neuronParamLowerBounds, neuronParamUpperBounds);

if exist(logFilePath, 'file') ~= 2
if exist(conciseFilePath, 'file') ~= 2

%}

