function m3ha_log_errors_params (outparams, err, simplexOut)
%% Log errors and parameter values
% Usage: m3ha_log_errors_params (outparams, err, simplexOut)
%
% Requires:
%       cd/save_params.m
%
% Used by:
%       ~/m3ha/optimizer4gabab/optimizer_4compgabab.m
%       ~/m3ha/optimizer4gabab/fminsearch3_4compgabab.m

% File History:
% 2017-01-17 Created
% 2017-01-25 Added concise version
% 2017-04-21 Added outparams.fitreg and outparams.fitregCpr
% 2017-05-15 Added currentparams to save to .p file
% 2017-05-16 Replaced tolx and tolf with relativeParamTolerance and relativeErrorTolerance
% 2017-05-17 Modified what to show and added outparams.logSwpErrFlag, 
%        outparams.logSwpWeightsFlag, outparams.logLtsWeightsFlag, 
%        outparams.logSwpEndsFlag
% 2017-05-17 If lts not computed, do not show errors
% 2018-08-10 Updated fitreg -> fitWindow
% TODO: Make each iteration a new column instead of a new row

% Fields to check if exist in err
tocheck_lts = {'avgltsverr', 'avgltsterr', 'avgltsdvdtverr'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if logging simplex version
if nargin < 3 || isempty(simplexOut)
    issimplex = 0;
else
    issimplex = 1;
end

% Extract info from outparams or simplexParams
numswps = outparams.numswps;
numparams = outparams.numparams;
if issimplex
    numinuse = length(simplexOut.paramsinuse);    
                    % number of parameters in use
    paramnamesinuse = simplexOut.paramnamesinuse;
    neuronparams = simplexOut.neuronparams;
    sheetName = fullfile(outparams.outFolderName, ...
        [outparams.prefix, ...
        '_simplexrun_', num2str(outparams.simplexnum), ...
        '_iter_', num2str(simplexOut.ctIterations), '_params.csv']);
else 
    neuronparams = outparams.neuronparams;
    sheetName = fullfile(outparams.outFolderName, ...
                        [outparams.prefix, '_params.csv']);
end
fitWindowCpr = outparams.fitWindowCpr;
fitWindow = outparams.fitWindow;
neuronparamnames = outparams.neuronparamnames;
neuronparams_min = outparams.neuronparams_min;
neuronparams_max = outparams.neuronparams_max;

% Check if LTS error exists
islts = isfield(err, tocheck_lts);

% Create file if not exist and print log header
fullfilename = fullfile(outparams.outFolderName, outparams.logfilename);
if exist(fullfilename, 'file') ~= 2
    fid = fopen(fullfilename, 'w');        % create log file for writing
    if issimplex
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
    if issimplex
        fprintf(fid, '%s, ', '# of parameters in use');
    end
    fprintf(fid, '%s, %s, ', ...
            'Total error', 'Normalized average sweep error');
    if all(islts)
        fprintf(fid, '%s, %s, %s, %s, ', ...
            'Normalized average LTS amp error', ...
            'Normalized average LTS time error', ...
            'Normalized average LTS slope error', ...
            'Normalized average LTS error');
    end
    fprintf(fid, repmat('%s, ', 1, numparams), neuronparamnames{:});
    if outparams.logSwpErrFlag
        strings = repmat({'sweep error (mV)'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
    end
    if all(islts) && outparams.logSwpErrFlag
        strings = repmat({'LTS amp error'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
        strings = repmat({'LTS time error'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
        strings = repmat({'LTS slope error'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
    end
    if outparams.logSwpWeightsFlag
        strings = repmat({'sweep weight'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
    end
    if all(islts) && outparams.logLtsWeightsFlag
        fprintf(fid, repmat('%s, ', 1, 3), ...
            'LTS amp weight', 'LTS time weight', ...
            'LTS slope weight');
        fprintf(fid, '%s, ', 'LTS to sweep error ratio');
    end
    fprintf(fid, '%s, %s, %s, %s, ', ...
        'cpr start', 'cpr end', 'sweep start', 'sweep end');
    %{
    if outparams.logSwpEndsFlag
        strings = repmat({'cpr start'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
        strings = repmat({'cpr end'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
        strings = repmat({'sweep start'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
        strings = repmat({'sweep end'}, 1, numswps);
        fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
    end
    %}
    fprintf(fid, '\n');
    fclose(fid);
end

% Log errors and params
fid = fopen(fullfilename, 'a');            % append
if issimplex
    fprintf(fid, '%6.4g, %6.4g, %s, %6.4g, %6.4g, %6.4g, ', ...
            simplexOut.ctIterations, simplexOut.maxIterations, ...
            simplexOut.lasthow, simplexOut.ctEvals, ...
            simplexOut.maxFunctionEvaluations, simplexOut.toterr);
    fprintf(fid, '%6.4g, %6.4g, %6.4g, %6.4g, ', ...
            simplexOut.maxErrorChange, simplexOut.relativeErrorTolerance, ...
            simplexOut.maxParamChange, simplexOut.relativeParamTolerance);
end
fprintf(fid, '%2.2g, ', numswps);
if issimplex
    fprintf(fid, '%2.2g, ', numinuse);
end
fprintf(fid, '%6.4g, %6.4g, ', err.toterr, err.avgswperr);
if all(islts)
    fprintf(fid, '%6.4g, %6.4g, %6.4g, %6.4g, ', ...
            err.avgltsverr, err.avgltsterr, ...
            err.avgltsdvdtverr, err.avgltserr);
end
fprintf(fid, repmat('%6.4g, ', 1, numparams), neuronparams);
if outparams.logSwpErrFlag
    fprintf(fid, repmat('%6.4g, ', 1, numswps), err.swperr);
end
if all(islts) && outparams.logSwpErrFlag
    fprintf(fid, repmat('%6.4g, ', 1, numswps), err.ltsverr);
    fprintf(fid, repmat('%6.4g, ', 1, numswps), err.ltsterr);
    fprintf(fid, repmat('%6.4g, ', 1, numswps), err.ltsdvdtverr);
end
if outparams.logSwpWeightsFlag
    fprintf(fid, repmat('%6.4g, ', 1, numswps), outparams.sweepWeights);
end
if all(islts) && outparams.logLtsWeightsFlag
    fprintf(fid, repmat('%6.4g, ', 1, 3), outparams.ltsWeights);
    fprintf(fid, '%6.4g, ', outparams.lts2SweepErrorRatio);
end
fprintf(fid, '%6.4g, ', fitWindowCpr(1));
fprintf(fid, '%6.4g, ', fitWindowCpr(2));
fprintf(fid, '%6.4g, ', fitWindow(1));
fprintf(fid, '%6.4g, ', fitWindow(2));
%{
if outparams.logSwpEndsFlag
    fprintf(fid, repmat('%6.4g, ', 1, numswps), fitWindowCpr(:, 1)');
    fprintf(fid, repmat('%6.4g, ', 1, numswps), fitWindowCpr(:, 2)');
    fprintf(fid, repmat('%6.4g, ', 1, numswps), fitWindow(:, 1)');
    fprintf(fid, repmat('%6.4g, ', 1, numswps), fitWindow(:, 2)');
end
%}
fprintf(fid, '\n');
fclose(fid);

% Create concise version if not exist and print log header
fullfilename_c = fullfile(outparams.outFolderName, ...
            outparams.logconcisefilename);
if exist(fullfilename_c, 'file') ~= 2
    fid = fopen(fullfilename_c, 'w');        % create
    if issimplex
        fprintf(fid, '%s, %s, %s, %s, ', ...
            'ctIterations', 'how', 'ctEvals', 'Total error');
        fprintf(fid, '%s, %s, ', ...
            'Maximum error change', 'Maximum parameter change');
    else
        fprintf(fid, '%s, ', 'Total error');
    end
    fprintf(fid, '%s, ', 'Normalized average sweep error');
    if all(islts)
        fprintf(fid, '%s, %s, %s, %s, ', ...
            'Normalized average LTS amp error', ...
            'Normalized average LTS time error', ...
            'Normalized average LTS slope error', ...
            'Normalized average LTS error');
    end
    if issimplex
        fprintf(fid, repmat('%s, ', 1, numinuse), paramnamesinuse{:});
    else
        fprintf(fid, repmat('%s, ', 1, numparams), neuronparamnames{:});
    end
    if all(islts) && outparams.logLtsWeightsFlag
        fprintf(fid, repmat('%s, ', 1, 3), ...
            'LTS amp weight', 'LTS time weight', ...
            'LTS slope weight');
        fprintf(fid, '%s, ', 'LTS to sweep error ratio');
    end
    fprintf(fid, '\n');
    fclose(fid);
end

% Log errors and params (concise version)
fid = fopen(fullfilename_c, 'a');        % append
if issimplex
    fprintf(fid, '%6.4g, %s, %6.4g, %6.4g, ', ...
            simplexOut.ctIterations, simplexOut.lasthow, ...
            simplexOut.ctEvals, simplexOut.toterr);
    fprintf(fid, '%6.4g, %6.4g, ', ...
            simplexOut.maxErrorChange, simplexOut.maxParamChange);
else
    fprintf(fid, '%6.4g, ', err.toterr);
end
fprintf(fid, '%6.4g, ', err.avgswperr);
if all(islts)
    fprintf(fid, '%6.4g, %6.4g, %6.4g, %6.4g, ', ...
            err.avgltsverr, err.avgltsterr, ...
            err.avgltsdvdtverr, err.avgltserr);
end
if issimplex
    fprintf(fid, repmat('%6.4g, ', 1, numinuse), simplexOut.paramsinuse);
else
    fprintf(fid, repmat('%6.4g, ', 1, numparams), neuronparams);
end
if all(islts) && outparams.logLtsWeightsFlag
    fprintf(fid, repmat('%6.4g, ', 1, 3), outparams.ltsWeights);
    fprintf(fid, '%6.4g, ', outparams.lts2SweepErrorRatio);
end
fprintf(fid, '\n');
fclose(fid);

%% Save parameters into a p file
%   Note: Always do this if not in simplex mode
if ~issimplex || outparams.saveParamsFlag
    save_params(neuronParamsTable, 'FileName', sheetName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE

fprintf(fid, '%2.2g, ', numswps);
fprintf(fid, '%6.4f, ', outparams.lts2SweepErrorRatio);
fprintf(fid, repmat('%6.4e, ', 1, numel(outparams.neuronparams)), outparams.neuronparams);
fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.swpw)), outparams.swpw);
fprintf(fid, '%6.4e, ', err.toterr);
fprintf(fid, repmat('%6.4e, ', 1, numel(err.swperr)), err.swperr);
fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.ltsWeights)), outparams.ltsWeights);
fprintf(fid, repmat('%6.4e, ', 1, 3), err.ltserrv, err.ltserrt, err.ltserrdvdtv);
fprintf(fid,'\n');

strings = repmat({'sweep left edge'}, 1, numswps);
fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
strings = repmat({'sweep right edge'}, 1, numswps);
fprintf(fid, repmat('%s, ', 1, numswps), strings{:});

fprintf(fid, repmat('%6.4g, ', 1, numswps), outparams.swpedges(:,1)');
fprintf(fid, repmat('%6.4g, ', 1, numswps), outparams.swpedges(:,2)');

fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.fitreg(:, 1))), outparams.fitreg(:, 1)');
fprintf(fid, repmat('%6.4f, ', 1, numel(outparams.fitreg(:, 2))), outparams.fitreg(:, 2)');

if outparams.cprflag

if issimplex
    fprintf(fid, repmat('%6.4g, ', 1, numparams), simplexOut.neuronparams);
else
    fprintf(fid, repmat('%6.4g, ', 1, numparams), outparams.neuronparams);
end

tocheck = {'ltserr_maxamp_val', 'ltserr_maxamp_time', 'ltserr_maxdiffamp_val'};
fprintf(fid, '%s, %s, %s', 'LTS amp error', 'LTS time error', 'LTS slope error');
fprintf(fid, '%6.4g, %6.4g, %6.4g, ', err.ltserr_maxamp_val, err.ltserr_maxamp_time, err.ltserr_maxdiffamp_val);

strings = repmat({'sweep error'}, 1, numswps);
fprintf(fid, repmat('%s, ', 1, numswps), strings{:});
fprintf(fid, repmat('%6.4g, ', 1, numswps), err.swperr);

for k = 1:numel(tocheck)
    if ~isfield(err, tocheck_lts{k})
        err.(tocheck_lts{k}) = NaN;
    end
end

islts = zeros(size(tocheck_lts));
for k = 1:numel(tocheck_lts)
    if isfield(err, tocheck_lts{k})
        islts(k) = 1;
    end
end

pfilename = fullfile(outparams.outFolderName, [outparams.prefix, '.p']);
%% Save parameters into a p file
currentparams = [neuronparamnames', num2cell(neuronparams'), ...
    num2cell(neuronparams_min'), num2cell(neuronparams_max')];
save(pfilename, 'currentparams');

pfilename = fullfile(outparams.outFolderName, ...
    [outparams.prefix, ...
    '_simplexrun_', num2str(outparams.simplexnum), ...
    '_iter_', num2str(simplexOut.ctIterations), '.p']);

currentparams = [neuronparamnames', num2cell(neuronparams'), ...
    num2cell(neuronparams_min'), num2cell(neuronparams_max')];
save(pfilename, 'currentparams');

save_params(sheetName, neuronparamnames, neuronparams, ...
            neuronparams_min, neuronparams_max);

%}

