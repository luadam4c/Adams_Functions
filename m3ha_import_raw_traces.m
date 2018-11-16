function [dCpr, dIpscr, thisInd, dCprAll] = m3ha_import_raw_traces (dataDir, nrows, nprow, g200pInd, pInd, rowInd, fnrow, vrow, cpstartwin, cprwinOrig, currentPulseAmplitudeOrig, mvw, actVhold, ipscTimeOrig, ipscDur, initialSlopes)
%% Import raw traces
%
% Requires:
%       cd/compute_rms_error.m
%       cd/correct_unbalanced_bridge.m
%       cd/find_ind_str_in_cell.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/singleneuronfitting42.m and later versions
%
% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-12-21 Changed tabs to spaces
% 2018-05-16 Now passes in vrow
% 2018-05-18 Changed all variables to camelback case
% 2018-05-18 Averaged the current pulse responses according to vHold
% 2018-05-21 Added outparams.baseNoiseIpscr and outparams.baseNoiseCpr
% 2018-06-20 Now uses initialSlopes to take 
%               out cpr traces with out-of-balance bridges
% 2018-06-21 Fixed bugs
% 2018-07-09 Added numswpsCpr as output
% 2018-07-09 Now uses compute_rms_error() instead of rms_Gaussian()
% 2018-07-31 Use correct_unbalanced_bridge.m to correct for out-of-balance bridges
% 2018-08-09 Now computes sweep weights here
% 2018-08-10 baseNoiseIpscr and baseNoiseCpr are now column vectors
% 2018-09-12 Added outparams.oldAverageCprFlag
% 2018-11-15 Moved to Adams_Functions

% Hard-coded constants
PA_PER_NA = 1000;

global outparams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions across servers
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for compute_rms_error.m, correct_unbalanced_bridge, 
    %   find_ind_str_in_cell.m, print_cellstr.m
    addpath(fullfile(functionsDirectory, 'Adams_Functions')); 
end

%% Prepare
% Find total number of sweeps to import
numswps = nrows * nprow;                        % total number of traces used

% Find the expected current pulse start time
cpstartExpected = mean(cpstartwin);

% Read from initialStopes
initialSlopeFilenames = initialSlopes.filenamesSorted;
initialSlopeThreshold1IndexBalanced = initialSlopes.iThreshold1Balanced;
initialSlopeThreshold2IndexBalanced = initialSlopes.iThreshold2Balanced;

% Create log file name
logFileName = fullfile(outparams.outFolderName, sprintf('%s.log', mfilename));
fid = fopen(logFileName, 'w');

%% Import sweeps
% Initialize output variables
dCpr = cell(numswps, 1);                        % stores current pulse response traces for fitting
dIpscr = cell(numswps, 1);                      % stores IPSC response traces for fitting
baseNoiseIpscr = zeros(numswps, 1);             % stores baseline noise in voltage recordings
filenames = cell(numswps, 1);                   % stores file names of traces used
currentPulseAmplitude = zeros(numswps, 1);      % stores the current pulse amplitude (nA)
% outparams.cpstart = zeros(numswps, 1);        % stores the time of current pulse application (ms)
holdPotentialIpscr = zeros(numswps, 1);         % stores the holding potentials right before IPSC application (mV)
ct = 0;                                         % counts number of raw traces imported
for rowi = 1:nrows
    for k = 1:nprow
        ct = ct + 1;
        if ~isempty(g200pInd{1})
            thisInd(ct) = g200pInd{rowi}(k);
        elseif ~isempty(pInd{1})
            thisInd(ct) = pInd{rowi}(k);
        elseif ~isempty(rowInd{1})
            thisInd(ct) = rowInd{rowi}(k);
        end
        filenames{ct} = fnrow{thisInd(ct)};
        m = matfile(fullfile(dataDir, '/matfiles/', filenames{ct}));
        fprintf('Using trace %s\n', filenames{ct});

        % Use original data
        tvecOrig = m.d_orig(:, 1);

        % Find current pulse response window
        %   Must be sure to include current pulse start and cpstartwinLength adjustment
        acprwinBegin = find(tvecOrig >= cprwinOrig(1), 1);
        cpstartwinLength = cpstartwin(2) - cpstartwin(1);
        acprwinEnd = find(tvecOrig <= cprwinOrig(2) + cpstartwinLength, 1, 'last');
        acprwinInd = acprwinBegin:acprwinEnd;       % indices corresponding to approximate cpr window

        % Use original data for current pulse response
        simsCpr = tvecOrig(2) - tvecOrig(1);        % sampling interval in ms
        tvecCpr = tvecOrig(acprwinInd);             % time vector of original data in ms
        gvecCpr = m.d_orig(acprwinInd, 2);          % conductance vector of original data in nS
        gvecCpr = gvecCpr / 1000;                   % conductance vector in uS
        ivecCpr = m.d_orig(acprwinInd, 3);          % current vector of original data in pA
        ivecCpr = ivecCpr / 1000;                   % current vector in nA
        vvecCpr = m.d_orig(acprwinInd, 4);          % voltage vector of original data in mV

        % Find the expected baseline window length in samples
        cprbasewinLength = round((cpstartExpected - cprwinOrig(1))/simsCpr);

        % Find current pulse amplitude (convert to nA) and start of current pulse application
        currentPulseAmplitude(ct, 1) = currentPulseAmplitudeOrig(thisInd(ct)) / PA_PER_NA;
        cpstartwinEnd = find(tvecCpr <= cpstartwin(2), 1, 'last');
        cpstartInd = find(ivecCpr(1:cpstartwinEnd) > currentPulseAmplitude(ct, 1) * 0.25, 1, 'last');
%         outparams.cpstart(ct, 1) = tvecCpr(cpstartInd);                 % current pulse start in ms

        % Find the number of indices to pad before current pulse response data
        cprstartInd = cpstartInd - cprbasewinLength;
        if cprstartInd < 1
            nIndToPadCprBase = 1 - cprstartInd;
        else
            nIndToPadCprBase = 0;
        end

        % Find indices for current pulse response relative to tvecCpr
        cprwinLength = round((cprwinOrig(2) - cprwinOrig(1))/simsCpr);
        cprendInd = cprstartInd + cprwinLength - 1;
        if cprstartInd < 1
            cprInd = 1:cprendInd;
        else
            cprInd = cprstartInd:cprendInd;
        end

        % Get full time vector of current pulse response
        tvecCprLength = round(outparams.cprWindow(2)/simsCpr);
        tvecCprFull = simsCpr * (1:tvecCprLength)';
        dCpr{ct}(:, 1) = tvecCprFull;                           % time vector of current pulse response

        % Pad raw traces so that current pulses lie at 2100-2110 ms
        nIndToPadBeforeCpr = round(outparams.cprWindow(1)/simsCpr);
        nIndToPad = nIndToPadCprBase + nIndToPadBeforeCpr;
        vvectopad = NaN * ones(1, nIndToPad)';
        ivectopad = zeros(1, nIndToPad)';
        gvectopad = zeros(1, nIndToPad)';
        dCpr{ct}(:, 2) = vertcat(vvectopad, vvecCpr(cprInd));   % voltage vector of current pulse response
        dCpr{ct}(:, 3) = vertcat(ivectopad, ivecCpr(cprInd));   % current vector of current pulse response
        dCpr{ct}(:, 4) = vertcat(gvectopad, gvecCpr(cprInd));   % conductance vector of current pulse response

        % Use median-filtered & resampled data for IPSC response
        tvec = m.d_mfrs(:, 1);  % time vector of median-filtered then resampled data in ms
        gvec = m.d_mfrs(:, 2);  % conductance vector of median-filtered then resampled data in nS
        gvec = gvec / 1000;     % conductance vector in uS
        ivec = m.d_mfrs(:, 3);  % current vector of median-filtered then resampled data in pA
        ivec = ivec / 1000;     % current vector in nA
        vvec = m.d_mfrs(:, 4);  % voltage vector of median-filtered then resampled data in mV

        % Holding potential was already extracted during data analysis
        holdPotentialIpscr(ct, 1) = actVhold(thisInd(ct));    

        % Pad raw traces so that IPSCs begin at 3000 ms (outparams.timeIPSC) 
        %            and goes on for 7000 ms (ipscDur)
        sims = tvec(2) - tvec(1);       % Should be 1 ms
        about3000 = round(outparams.timeIPSC/sims)*sims;
        tvectopad = (sims:sims:about3000)';
        vvectopad = holdPotentialIpscr(ct, 1) * ones(1, length(tvectopad))';
        ivectopad = zeros(1, length(tvectopad))';
        gvectopad = zeros(1, length(tvectopad))';
        indofipsc = round(ipscTimeOrig/sims);
        indofend = round((ipscTimeOrig + ipscDur)/sims);
        tvec_shifted = tvec(indofipsc:indofend) - tvec(indofipsc - 1) + about3000;
        vvec_shifted = vvec(indofipsc:indofend);
        ivec_shifted = ivec(indofipsc:indofend);
        gvec_shifted = gvec(indofipsc:indofend);
        dIpscr{ct}(:, 1) = vertcat(tvectopad, tvec_shifted);        % time vector of IPSC response
        dIpscr{ct}(:, 2) = vertcat(vvectopad, vvec_shifted);        % voltage vector of IPSC response
        dIpscr{ct}(:, 3) = vertcat(ivectopad, ivec_shifted);        % current vector of IPSC response
        dIpscr{ct}(:, 4) = vertcat(gvectopad, gvec_shifted);        % conductance vector of IPSC response

        % Find the baseline noise
        indBaseline = 1:indofipsc-1;
        baseline = vvec(indBaseline);
        baseNoiseIpscr(ct) = compute_rms_error(baseline);
    end
end
if ct ~= numswps
    error('Number of traces imported incorrect!!');
else
    fprintf('\n');
end

%% Fix current pulse response traces that may have out-of-balance bridges
if outparams.correctDcStepsFlag
    % Find the index of file in all files sorted by initial slope
    %   in descending order
    ftemp = @(x) find_ind_str_in_cell(x, initialSlopeFilenames);

    % Determine whether the initial slopes exceed threshold
    %   Note: These may have out-of-balance bridges
    isOutOfBalance = ...
        cellfun(@(x) ftemp(x) < initialSlopeThreshold2IndexBalanced || ...
                    ftemp(x) > initialSlopeThreshold1IndexBalanced, ...
                    filenames);

    % Print out an appropriate message
    if any(isOutOfBalance)
        fprintf(fid, ['The following current pulse responses will be ', ...
                        'corrected due to out-of-balance bridges:\n']);
        print_cellstr(filenames(isOutOfBalance), ...
                      'FileID', fid, 'OmitBraces', true, ...
                      'OmitQuotes', true, 'Delimiter', '\n');
    else
        fprintf(fid, ['There are no current pulse response traces ', ...
                        'with out-of-balance bridges!\n']);
    end

    % Correct for traces that may have out-of-balance bridges
    parfor iSwp = 1:numswps
        if isOutOfBalance(iSwp)
            % Get the old trace
            vvecOld = dCpr{iSwp}(:, 2);
            ivecOld = dCpr{iSwp}(:, 3);

            % Correct any unbalanced bridge in the trace
            vvecNew = correct_unbalanced_bridge(vvecOld, ivecOld);

            % Store the new trace
            dCpr{iSwp}(:, 2) = vvecNew;
        end
    end
end

%% Average the current pulse responses according to vHold
if outparams.oldAverageCprFlag
    % Find the holding voltage assigned for each sweep
    vholdThis = vrow(thisInd);

    % Get the number of data points for the current pulse response
    ndpCpr = size(dCpr{1}, 1);

    % Find unique vHold values
    vUnique = sort(unique(vholdThis));
    nVhold = length(vUnique);

    % Group the traces by unique vHold values, then average the grouped traces
    dataCprGrouped = cell(nVhold, 1);
    dataCprAveraged = cell(nVhold, 1);
    for iVhold = 1:nVhold
        % Get the current vHold value
        vnow = vUnique(iVhold);

        % Collect all cpr traces with this vHold value
        dCprGroupedThisVhold = dCpr(vholdThis == vnow);

        % Preallocate
        dCprAveragedThisVhold = zeros(ndpCpr, 4);

        % Take the time, current and conductance traces from the first trace
        dCprAveragedThisVhold(:, 1) = dCprGroupedThisVhold{1}(:, 1);
        dCprAveragedThisVhold(:, 3) = dCprGroupedThisVhold{1}(:, 3);
        dCprAveragedThisVhold(:, 4) = dCprGroupedThisVhold{1}(:, 4);

        % Average the voltage traces
        temp1 = cellfun(@(x) x(:, 2), dCprGroupedThisVhold, ...
                        'UniformOutput', false);
        vCprGroupedThisVhold = cell2mat(temp1');
        vCprAveragedThis = nanmean(vCprGroupedThisVhold, 2);
        dCprAveragedThisVhold(:, 2) = vCprAveragedThis;

        % Save in arrays
        dataCprGrouped{iVhold} = dCprGroupedThisVhold;
        dataCprAveraged{iVhold} = dCprAveragedThisVhold;
    end
else
    dataCprAveraged = dCpr;
end

%% Find the holding potential, holding current and baseline noise
% Count the number of sweeps for the current pulse response
numswpsCpr = numel(dataCprAveraged);

% Find the holding potential, holding current and baseline noise
baseNoiseCpr = zeros(numswpsCpr, 1);
holdCurrentCpr = zeros(numswpsCpr, 1);
holdPotentialCpr = zeros(numswpsCpr, 1);
for iSwp = 1:numswpsCpr
    % Get the current voltage trace
    vCprAveragedThis = dataCprAveraged{iSwp}(:, 2);

    % Find the baseline noise
    indBaselineCpr = nIndToPadBeforeCpr + (1:cprbasewinLength);
    baselineCprThis = vCprAveragedThis(indBaselineCpr);
    cprBaselineNoise = compute_rms_error(baselineCprThis);

    % Find holding potential right before current pulse application
    mvind = 1:round(mvw/simsCpr);   % indices for taking the mean of voltages
    baseInd = nIndToPadBeforeCpr + cpstartInd - ...
                round(0.5/simsCpr) - fliplr(mvind);      % base indices
    cprBaselineVoltage = nanmean(vCprAveragedThis(baseInd));

    % Estimate the holding currents (nA) based on estimated
    %   input resistance and resting membrane potential
    %   I = (V - epas) / R
    cprHoldCurrent = (cprBaselineVoltage - epasEstimate) / RinEstimate;

    % Save in arrays
    holdPotentialCpr(iSwp) = cprBaselineVoltage;
    holdCurrentCpr(iSwp) = cprHoldCurrent;
    baseNoiseCpr(iSwp) = cprBaselineNoise;
end

if outparams.generateDataFlag
    % Estimate the corresponding variations in holding current (nA)
    holdCurrentNoiseIpscr = (baseNoiseIpscr / RinEstimate) * 5;
    holdCurrentNoiseCpr = (baseNoiseCpr / RinEstimate) * 5;
else
    holdCurrentNoiseIpscr = (baseNoiseIpscr / RinEstimate) * 0;
    holdCurrentNoiseCpr = (baseNoiseCpr / RinEstimate) * 0;
end

% Compute sweep weights based on baseline noise
sweepWeightsCpr = 1 ./ baseNoiseCpr;
sweepWeightsIpscr = 1 ./ baseNoiseIpscr;

%% Wrap up
% Rename as new dCpr
dCprAll = dCpr;
dCpr = dataCprAveraged;

% Save in outparams
outparams.filenames = filenames;
outparams.currentPulseAmplitude = currentPulseAmplitude;

outparams.holdPotentialCpr = holdPotentialCpr;
outparams.holdCurrentCpr = holdCurrentCpr;
outparams.baseNoiseCpr = baseNoiseCpr;
outparams.holdCurrentNoiseCpr = holdCurrentNoiseCpr;
outparams.sweepWeightsCpr = sweepWeightsCpr;

outparams.holdPotentialIpscr = holdPotentialIpscr;
outparams.holdCurrentIpscr = holdCurrentIpscr;
outparams.baseNoiseIpscr = baseNoiseIpscr;
outparams.holdCurrentNoiseIpscr = holdCurrentNoiseIpscr;
outparams.sweepWeightsIpscr = sweepWeightsIpscr;

% Close log file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find approximate current pulse response window
% Must be sure to include current pulse start and at least 500 ms after current pulse start
acprwinBegin = find(m.d_orig(:, 1) >= cpstartwin(1), 1);
cpstartwinLength = cpstartwin(2) - cpstartwin(1);
acprwinEnd = find(m.d_orig(:, 1) <= cprwinOrig(2) + cpstartwinLength/2, 1, 'last');
acprwinInd = acprwinBegin:acprwinEnd;       % indices corresponding to approximate cpr window

% Find indices for current pulse response relative to tvecCpr
cprstartInd = cpstartInd;
cprwinLength = cprwinOrig(2) - cprwinOrig(1);
cprendInd = cprstartInd + round(cprwinLength/simsCpr);
cprInd = (cprstartInd + 1):cprendInd;       % since we are padding later, we need the +1

% Find holding potential right before current pulse application
mvind = 1:round(mvw/simsCpr);               % indices for taking the mean of voltages
baseInd = cpstartInd - round(0.5/simsCpr) - fliplr(mvind);      % base indices
cprBaselineVoltage = mean(vvecCpr(baseInd));
outparams.holdPotentialCpr(ct, 1) = cprBaselineVoltage;

% Pad raw traces so that current pulses lie at 2100-2110 ms and response goes on until 2500 ms
about2100 = round(outparams.cprWindow(1)/simsCpr) * simsCpr;       % about 2100 ms
                % this is needed because simsCpr is often 0.099 ms
tvectopad = (simsCpr:simsCpr:about2100)';
vvectopad = outparams.holdPotentialCpr(ct, 1) * ones(1, length(tvectopad))';
ivectopad = zeros(1, length(tvectopad))';
gvectopad = zeros(1, length(tvectopad))';
tvec_shifted = tvecCpr(cprInd) - tvecCpr(cprstartInd - 1) + about2100;
vvec_shifted = vvecCpr(cprInd);
ivec_shifted = ivecCpr(cprInd);
gvec_shifted = gvecCpr(cprInd);
dCpr{ct}(:, 1) = vertcat(tvectopad, tvec_shifted);    % time vector of current pulse response
dCpr{ct}(:, 2) = vertcat(vvectopad, vvec_shifted);    % voltage vector of current pulse response
dCpr{ct}(:, 3) = vertcat(ivectopad, ivec_shifted);    % current vector of current pulse response
dCpr{ct}(:, 4) = vertcat(gvectopad, gvec_shifted);    % conductance vector of current pulse response

holdPotentialCpr = zeros(numswps, 1);          % stores the holding potentials right before cp application (mV)

%       /home/Matlab/Kojis_Functions/rms_Gaussian.m
    addpath(fullfile(functionsDirectory, 'Kojis_Functions')); 
                                            % for rms_Gaussian.m
    baseNoiseIpscr(ct) = rms_Gaussian(baseline);
    cprBaselineNoise = rms_Gaussian(baselineCprThis);

% Take out traces that may have out-of-balance bridges if requested
if outparams.correctDcStepsFlag
    % Determine whether the initial slopes exceed threshold
    %   Note: These may have out-of-balance bridges
    isOutOfBalance = ...
        cellfun(@(x) find_ind_str_in_cell(x, initialSlopeFilenames) < ...
                                        initialSlopeThresholdIndexToInclude, ...
                    fileNamesThisVhold);

    % Print out appropriate message
    if any(isOutOfBalance)
        fprintf(fid, ['For Vhold == %g, ', ...
                 'the following current pulse responses will be taken out ', ...
                 'due to out-of-balance bridges:\n'], vnow);
        print_cellstr(fileNamesThisVhold(isOutOfBalance), ...
                      'FileID', fid, 'OmitBraces', true, ...
                      'OmitQuotes', true, 'Delimiter', '\n');
    else
        fprintf(fid, ['For Vhold == %g, ', ...
                 'there are no out-of-balance bridges!\n'], vnow);
    end

    % Take out traces that may have out-of-balance bridges
    dCprGroupedThisVhold = dCprGroupedThisVhold(~isOutOfBalance);
end

initialSlopeThresholdIndexToInclude = initialSlopes.iThreshold1Balanced;

holdCurrentIpscr = zeros(numswps, 1);                % stores the holding currents recorded in Christine's Excel sheet (nA)

% Average the holding currents
holdCurrentThisVhold = mean(holdCurrentIpscr(vholdThis == vnow));

baseNoiseCpr = zeros(1, nVhold);

baseNoiseIpscr = zeros(1, numswps);                  % stores baseline noise in voltage recordings

isOutOfBalance = ...
    cellfun(@(x) ftemp(x) < initialSlopeThreshold2IndexBalanced || ...
                ftemp(x) > initialSlopeThreshold1IndexBalanced, ...
                fileNamesThisVhold);

% Correct for traces that may have out-of-balance bridges
parfor iTrace = 1:nTracesThisVhold
    if isOutOfBalance(iTrace)
        % Get the old trace
        vvecOld = dCprGroupedThisVhold{iTrace}(:, 2);
        ivecOld = dCprGroupedThisVhold{iTrace}(:, 3);

        % Correct any unbalanced bridge in the trace
        vvecNew = correct_unbalanced_bridge(vvecOld, ivecOld);

        % Store the new trace
        dCprGroupedThisVhold{iTrace}(:, 2) = vvecNew;
    end
end

print_cellstr(fileNamesThisVhold(isOutOfBalance), ...
              'FileID', fid, 'OmitBraces', true, ...
              'OmitQuotes', true, 'Delimiter', '\n');

% Get the corresponding filenames
fileNamesThisVhold = filenames(vholdThis == vnow);


% Number of current pulse response traces to fit
numswpsCpr = nVhold;

% Find the number of traces
nTracesThisVhold = numel(dCprGroupedThisVhold);

%% Get other information for each sweep
% Find the holding current used by Christine for each sweep, 
%   converting from pA to nA
holdCurrentIpscr = actIhold(thisInd)' / 1000;

% Compute the regression coefficient
% RinRegression = (holdCurrentIpscr - mean(holdCurrentIpscr)) \ ...
%                     (holdPotentialIpscr - mean(holdPotentialIpscr));

%% Estimate epas and Rin with linear least squares
% Define the matrices
% Ohm's Law: IR + epas = V
%     units: [nA] * [MOhm] + [mV] = [mV]
% X * w = V
X = ones(numswps, 2);
X(:, 1) = holdCurrentIpscr;
V = holdPotentialIpscr;

% Compute the estimates
w = pinv(X) * V;
RinEstimate = w(1);                 % input resistance (MOhm)
epasEstimate = w(2);                % resting membrane potential (mV)

% Check if the values make sense
if RinEstimate <= 0
    colorRin = 'r';
else
    colorRin = 'k';
end

% Construct a vector of holding currents
holdCurrentToPlot = linspace(min(holdCurrentIpscr), max(holdCurrentIpscr), 1000);

% Compute predicted values
holdPotentialPredicted = holdCurrentToPlot * RinEstimate + epasEstimate;

% Plot holding potential versus holding current
h = figure('Visible', 'off');
clf(h);
hold on;
plot(holdCurrentIpscr, holdPotentialIpscr, 'o', 'LineWidth', 2);
plot(holdCurrentToPlot, holdPotentialPredicted, 'r')
text(0.1, 0.9, ['Rin = ', num2str(RinEstimate), ' MOhm'], ...
    'Units', 'normalized', 'Color', colorRin);
text(0.1, 0.85, ['epas = ', num2str(epasEstimate), ' mV'], ...
    'Units', 'normalized', 'Color', 'k');
% text(0.1, 0.8, ['Slope = ', num2str(RinRegression), ' MOhm'], ...
%     'Units', 'normalized', 'Color', 'r');
title(['Voltage-Current relationship for ', outparams.cellName]);
ylabel('Holding potential (mV)');
xlabel('Holding current (nA)');
figName = fullfile(outparams.outFolderName, ...
                    [outparams.prefix, ...
                    '_voltage-current-relationship.png']);
saveas(h, figName);
close(h)


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%