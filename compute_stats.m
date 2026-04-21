function varargout = compute_stats (vecs, statName, varargin)
%% Computes a statistic of vector(s) possibly restricted by endpoint(s)
% Usage: [stat1, stat2, ...] = compute_stats (vecs, statName, dim (opt), varargin)
% Explanation:
%       Computes the specified statistic(s) for the provided vectors or cell array 
%       of vectors. Supported statistics include measures of central tendency 
%       (mean, median), dispersion (std, stderr, quartiles, range), and 
%       intervals (confidence intervals). Can optionally ignore NaNs, remove 
%       outliers, or operate on specific sub-indices/windows.
%
%       Note: If any element is empty, returns NaN.
%
% Example(s):
%       data = randn(10, 3);
%       compute_stats(data, 'mean')
%       compute_stats(data, 'std')
%       [meanData, stdData] = compute_stats(data, {'mean', 'std'})
%       compute_stats(data, 'stderr')
%       compute_stats(data, 'err')
%       compute_stats(data, 'err95')
%       compute_stats(data, 'lower95')
%       compute_stats(data, 'upper95')
%       compute_stats(data, 'lower95med')
%       compute_stats(data, 'upper95med')
%       compute_stats(data, 'cov')
%       compute_stats(data, 'zscore')
%       compute_stats(data, 'mean', 2)
%       compute_stats(data, 'max', 2)
%       compute_stats(data, 'mean', 'IgnoreNan', true)
%
% Outputs:
%       varargout   - the computed statistic(s) for each vector
%                   specified as a numeric vector 
% Arguments:
%       vecs        - vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%       statName    - name(s) of the statistic
%                   must be a string/char or cell array of strings unambiguously matching: 
%                       'average' or 'mean' - mean
%                       'median'    - median
%                       'quartile25' - 25th percentile (first quartile)
%                       'quartile75' - 75th percentile (third quartile)
%                       'std'       - standard deviation
%                       'stderr'    - standard error
%                       'err' or 'err95' - error margin for the 95% confidence interval
%                       'lower95'   - lower bound of the 95% confidence interval of the mean
%                       'upper95'   - upper bound of the 95% confidence interval of the mean
%                       'lower95med' - lower bound of the 95% confidence interval of the median
%                       'upper95med' - upper bound of the 95% confidence interval of the median
%                       'cov'       - coefficient of variation
%                       'zscore'    - z-score
%                       'max'       - maximum
%                       'min'       - manimum
%                       'range'     - range
%                       'range2mean'- percentage of range relative to mean
%       dim         - (opt) dimension to compute stats along
%                   must be either 1, 2 or 3
%                   default == 1
%       varargin    - 'DataType': type of data to calculate stats for
%                   must be an unambiguous, case-insensitive match to one of:
%                       'normal'        - standard calculation
%                       'proportion'    - calculation for a proportion
%                       'ratio'         - calculation for a ratio
%                   default == 'normal'
%                   - 'Weights': optional weight vector/array for calculating 
%                       weighted statistics (mean, std, stderr, err, cov, zscore)
%                   must be empty, a numeric array, or a cell array of arrays
%                   default == []
%                   - 'IgnoreNan': whether to ignore NaN entries
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveOutliers': whether to remove outliers
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Indices': indices for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == set by extract_subvectors.m
%                   - 'Endpoints': endpoints for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == set by extract_subvectors.m
%                   - 'Windows': value windows to extract 
%                       Note: this assumes that the values are nondecreasing
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == set by extract_subvectors.m
%                   
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/isemptycell.m
%       cd/remove_outliers.m
%       cd/stderr.m
%       cd/nanstderr.m
%
% Used by:
%       cd/compute_combined_array.m
%       cd/compute_combined_trace.m
%       cd/compute_population_average.m
%       cd/compute_sampsizepwr.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/parse_multiunit.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/plot_autocorrelogram.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_chevron.m
%       cd/plot_chevron_bar_inset.m
%       cd/plot_grouped_jitter.m
%       cd/plot_traces.m
%       cd/select_similar_values.m
%       cd/test_difference.m
%       vIRt-Moore/jm_postprocess_virt_sim_ExtCurrent_analysis.m
%       scAAV/analyze_qupath_figure2.m
%       scAAV/analyze_qupath_figure3_spinal_cord.m
%
% Related functions:
%       cd/compute_weighted_average.m

% File History:
% 2018-12-17 Created by Adam Lu
% 2019-03-14 compute_means -> compute_stats
% 2019-03-14 Added statName as a required argument
% 2019-03-14 Added 'IgnoreNan' as an optional argument
% 2019-03-14 Added 'RemoveOutliers' as an optional argument
% 2019-03-14 Added 'cov' to validStatNames
% 2019-05-12 Added 'zscore', 'range' and 'range2mean'
% 2019-05-12 Added dim as an optional argument
% 2019-08-07 Fixed 0.95 -> 0.96
% 2019-08-20 Now always return NaN if empty
% 2019-09-19 Added 'max' and 'min'
% 2019-11-14 Fixed usage of std and nanstd
% 2019-11-26 Updated confidence intervals to use t-distribution
% 2019-11-27 Added 'err'
% 2026-01-09 Added 'median', 'quartile25', 'quartile75'
% 2026-01-09 Added 'lower95med' and 'upper95med' (CI of median)
% 2026-03-17 Added 'DataType' and 'Weights' optional arguments by Gemini
% 2026-03-17 Optimized nanstderr calls by Gemini
% 2026-03-17 Refactored logic into documented helper functions by Gemini
% 2026-03-17 Now optionally accept a cell array of statNames with varargout
% 2026-03-20 Implemented Wilson Score interval for proportion bounds by Gemini
% 2026-03-20 Fixed ratio datatype transformations, asymmetric bounds, and t-distribution degrees of freedom by Gemini
% TODO: Combine with compute_weighted_average.m

%% Hard-coded parameters
% Define the list of allowable statistic names to prevent invalid inputs
validStatNames = {'average', 'mean', 'median', ...
                    'quartile25', 'quartile75', ...
                    'std', 'stderr', 'err', 'err95', ...
                    'lower95', 'upper95', ...
                    'lower95med', 'upper95med', ...
                    'cov', 'zscore', 'max', 'min', 'range', 'range2mean'};

% Define the allowable data types for metric computations
validDataTypes = {'normal', 'proportion', 'ratio'};

%% Default values for optional arguments
dimDefault = 1;                 % compute across rows by default
dataTypeDefault = 'normal';     % assume standard numeric calculations by default
weightsDefault = [];            % no weights by default
ignoreNanDefault = false;       % don't ignore NaN by default
removeOutliersDefault = false;  % don't remove outliers by default
indicesDefault = [];            % extract entire trace(s) by default
endPointsDefault = [];          % extract entire trace(s) by default
windowsDefault = [];            % extract entire trace(s) by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vecs', ...                  % vectors to extract
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'statName', ...              % requested statistic
    @(x) ischar(x) || isstring(x) || iscellstr(x));

% Add optional inputs to the Input Parser
addOptional(iP, 'dim', dimDefault, ...       % dimension to compute along
    @(x) assert(isnumeric(x) && (x == 1 || x == 2 || x == 3), ...
                'dim must be either 1, 2 or 3!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataType', dataTypeDefault, ...
    @(x) any(validatestring(x, validDataTypes)));
addParameter(iP, 'Weights', weightsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                'Weights must be empty, a numeric array, or a cell array of numeric arrays!'));
addParameter(iP, 'IgnoreNan', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveOutliers', removeOutliersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Indices', indicesDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Indices must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Windows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'Windows', windowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Windows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vecs, statName, varargin{:});

% Extract standard arguments
dim = iP.Results.dim;
dataType = validatestring(iP.Results.DataType, validDataTypes);
weights = iP.Results.Weights;

% Extract logical flags
ignoreNan = iP.Results.IgnoreNan;
removeOutliers = iP.Results.RemoveOutliers;

% Extract subvector configurations
indices = iP.Results.Indices;
endPoints = iP.Results.EndPoints;
windows = iP.Results.Windows;

%% Preparation
% Convert statName to a cell array of strings if it is passed as a single string/char
if ischar(statName) || isstring(statName)
    statNames = {char(statName)};
else
    statNames = statName;
end

% Validate all requested statistic names against the hardcoded valid list
statNames = cellfun(@(x) validatestring(x, validStatNames), statNames, 'UniformOutput', false);

%% Do the job
% Extract specified subvectors if bounds or windows were passed
subVecs = extract_subvectors(vecs, 'Indices', indices, ...
                            'EndPoints', endPoints, 'Windows', windows);

% Extract subweights parallel to subvectors if weights are provided
if ~isempty(weights)
    % Outlier removal would cause weight misalignment, so prevent it
    if removeOutliers
        error('RemoveOutliers cannot be true when Weights are provided, as it may cause misalignment.');
    end
    % Apply identical subvector extraction logic to the weights
    subWeights = extract_subvectors(weights, 'Indices', indices, ...
                                'EndPoints', endPoints, 'Windows', windows);
else
    % Assign empty array if no weights exist
    subWeights = [];
end

% Remove outliers from the extracted subvectors if requested
if removeOutliers
    if iscell(subVecs)
        % Apply outlier removal to each cell independently
        subVecs = cellfun(@remove_outliers, subVecs, ...
                        'UniformOutput', false);
    else
        % Apply outlier removal to the array
        subVecs = remove_outliers(subVecs);
    end
end

% Initialize output cell array to dynamically hold multiple statistics
numStats = numel(statNames);
varargout = cell(1, numStats);

% Compute each requested statistic
for iStat = 1:numStats
    % Extract the name of the statistic currently being computed
    currentStatName = statNames{iStat};
    
    % Decide on the function to use on each vector for current statistic
    switch currentStatName
        case {'average', 'mean'}
            % Compute mean, handling data types, missing values, and weights
            func = @(x, w) compute_mean_helper(x, dim, dataType, ignoreNan, w);
        case 'median'
            % Compute median, enforcing unweighted fallback
            func = @(x, w) handle_unweighted_only(x, w, currentStatName, @(x) compute_median_helper(x, dim, ignoreNan));
        case {'quartile25', 'quartile75'}
            % Configure the correct percentile (25 or 75)
            p = 25; if strcmpi(currentStatName, 'quartile75'); p = 75; end
            
            % Compute quartile, enforcing unweighted fallback
            func = @(x, w) handle_unweighted_only(x, w, currentStatName, @(x) compute_quartile_helper(x, p, dim, ignoreNan));
        case 'std'
            % Compute standard deviation, passing dataType to correctly apply Delta Method for ratios
            func = @(x, w) compute_std_helper(x, dim, dataType, ignoreNan, w);
        case 'stderr'
            % Compute standard error, passing dataType to correctly apply Delta Method for ratios
            func = @(x, w) compute_stderr_helper(x, dim, dataType, ignoreNan, w);
        case {'err', 'err95', 'lower95', 'upper95'}
            % Define a helper to calculate the effective sample size for confidence intervals
            nFunc = @(x, w) compute_n_eff(x, w, ignoreNan, dim, dataType);
            
            % Set the appropriate critical value function (z for proportion, t otherwise)
            if strcmpi(dataType, 'proportion')
                % Use standard normal z-score for proportion distributions (1.95996)
                critFunc = @(x, w) norminv(0.975) .* ones(size(nFunc(x, w)));
            else
                % Use t-distribution for continuous data (approximating df as n-1)
                critFunc = @(x, w) arrayfun(@(y) tinv(0.975, max(1, y - 1)), nFunc(x, w));
            end

            % Define functions for computing boundaries and margins based on the specific stat request
            switch currentStatName
                case {'err', 'err95'}
                    % Use half-width for asymmetric distributions like proportion and ratio
                    if strcmpi(dataType, 'proportion') || strcmpi(dataType, 'ratio')
                        func = @(x, w) (compute_ci_bound_helper(x, dim, dataType, critFunc(x, w), 'upper', w, ignoreNan) - ...
                                        compute_ci_bound_helper(x, dim, dataType, critFunc(x, w), 'lower', w, ignoreNan)) / 2;
                    else
                        % Use symmetric multiplication for standard normal distributions
                        func = @(x, w) critFunc(x, w) .* compute_stderr_helper(x, dim, dataType, ignoreNan, w);
                    end
                case 'lower95'
                    % Compute the lower bound explicitly
                    func = @(x, w) compute_ci_bound_helper(x, dim, dataType, critFunc(x, w), 'lower', w, ignoreNan);
                case 'upper95'
                    % Compute the upper bound explicitly
                    func = @(x, w) compute_ci_bound_helper(x, dim, dataType, critFunc(x, w), 'upper', w, ignoreNan);
            end
        case {'lower95med', 'upper95med'}
            % Compute median bounds using an unweighted asymptotic approximation
            func = @(x, w) handle_unweighted_only(x, w, currentStatName, @(x) compute_ci_med_helper(x, dim, ignoreNan, currentStatName));
        case 'cov'
            % Compute coefficient of variation (std / mean), passing dataType for proper ratio handling
            func = @(x, w) compute_std_helper(x, dim, dataType, ignoreNan, w) ./ compute_mean_helper(x, dim, dataType, ignoreNan, w);
        case 'zscore'
            % Compute z-score (mean / std), passing dataType for proper ratio handling
            func = @(x, w) compute_mean_helper(x, dim, dataType, ignoreNan, w) ./ compute_std_helper(x, dim, dataType, ignoreNan, w);
        case 'max'
            % Compute maximum, enforcing unweighted fallback
            func = @(x, w) handle_unweighted_only(x, w, currentStatName, @(x) max(x, [], dim));
        case 'min'
            % Compute minimum, enforcing unweighted fallback
            func = @(x, w) handle_unweighted_only(x, w, currentStatName, @(x) min(x, [], dim));
        case 'range'
            % Compute range, enforcing unweighted fallback
            func = @(x, w) handle_unweighted_only(x, w, currentStatName, @(x) range(x, dim));
        case 'range2mean'
            % Compute range as a percentage of the mean, enforcing unweighted fallback for the range
            func = @(x, w) handle_unweighted_only(x, w, currentStatName, @(x) (range(x, dim) ./ compute_mean_helper(x, dim, dataType, ignoreNan, [])) * 100);
        otherwise
            error('Code logic error!');
    end

    % Process the computed statistic across the structure of subVecs
    if iscell(subVecs)
        % Check if the cell array is free of empty arrays
        if ~any(isemptycell(subVecs))
            if isempty(subWeights)
                % Apply unweighted function across all cells
                currentStats = cellfun(@(x) func(x, []), subVecs);
            else
                % Apply weighted function across paired data and weight cells
                currentStats = cellfun(func, subVecs, subWeights);
            end
        else
            % Pre-allocate a NaN array parallel to the cell array
            currentStats = nan(size(subVecs));
            
            % Iterate manually to safely evaluate non-empty cells
            parfor iVec = 1:numel(subVecs)
                if ~isempty(subVecs{iVec})
                    if isempty(subWeights)
                        % Compute stat for unweighted data
                        currentStats(iVec) = func(subVecs{iVec}, []);
                    else
                        % Compute stat for weighted data
                        currentStats(iVec) = func(subVecs{iVec}, subWeights{iVec});
                    end
                end
            end
        end
    else
        % For numeric arrays, evaluate directly or return NaN if empty
        if isempty(subVecs)
            currentStats = NaN;
        else
            currentStats = func(subVecs, subWeights);
        end
    end
    
    % Assign the finalized statistic array to the varargout cell
    varargout{iStat} = currentStats;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = compute_prctile_strict(x, p, dim)
% Helper to compute percentiles while propagating NaNs (strict behavior)

% Calculate standard percentiles using built-in function
y = prctile(x, p, dim);

% Identify any slices that contain at least one NaN
nans = any(isnan(x), dim);

% Overwrite percentile results with NaN where missing values were present
y(nans) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function statVal = get_nanstderr_stat_val(x, dim, dataType, w)
% Helper to efficiently extract the calculated stat value from nanstderr

% Assign empty weights if none are provided
if nargin < 4
    w = []; 
end

% Request both standard error and stat value, but only keep the stat value
[~, statVal] = nanstderr(x, dim, 'DataType', dataType, 'Weights', w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bound = compute_ci_bound_helper(x, dim, dataType, critVal, direction, w, ignoreNan)
% Helper to compute confidence bounds consolidating weighted logic

% Assign empty weights if none are provided
if nargin < 6
    w = []; 
end

% Process ratio data in the log domain and exponentiate to ensure strictly positive bounds
if strcmpi(dataType, 'ratio')
    % Log transform for computation
    xLog = log(x);
    
    % Set values that were 0 or negative to NaN since log is undefined there
    xLog(x <= 0) = NaN;
    
    % Compute the normal mean and standard error in the log space
    meanLog = compute_mean_helper(xLog, dim, 'normal', ignoreNan, w);
    stderrLog = compute_stderr_helper(xLog, dim, 'normal', ignoreNan, w);
    
    % Calculate the interval bounds in log space, then exponentiate back
    if strcmpi(direction, 'lower')
        bound = exp(meanLog - critVal .* stderrLog);
    else
        bound = exp(meanLog + critVal .* stderrLog);
    end
    
    % Apply standard NaN propagation if not strictly ignored
    if ~ignoreNan
        bound(any(isnan(x), dim)) = NaN;
    end
    
    % Return early
    return;
end

% Route to Wilson Score interval if the data type is a proportion to keep bounds between 0 and 1
if strcmpi(dataType, 'proportion')
    bound = compute_wilson_score_helper(x, dim, w, ignoreNan, critVal, direction);
    return;
end

% Check if we need to handle special data types, missing values, or weights
if ~strcmpi(dataType, 'normal') || ignoreNan || ~isempty(w)
    % Compute both standard error and the central statistic using nanstderr
    [stdErr, statVal] = nanstderr(x, dim, 'DataType', dataType, 'Weights', w);
    
    % Apply the critical value margin to the correct direction
    if strcmpi(direction, 'lower')
        bound = statVal - critVal .* stdErr;
    else
        bound = statVal + critVal .* stdErr;
    end
else
    % For standard, unweighted data with no missing values, use fast built-in functions
    if strcmpi(direction, 'lower')
        bound = mean(x, dim) - critVal .* stderr(x, dim);
    else
        bound = mean(x, dim) + critVal .* stderr(x, dim);
    end
end

% Propagate NaNs if not explicitly ignored
if ~ignoreNan
    bound(any(isnan(x), dim)) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bound = compute_wilson_score_helper(x, dim, w, ignoreNan, zVal, direction)
% Helper to compute the Wilson Score interval for proportions
% Reference: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

% Calculate the effective sample size
n = compute_n_eff(x, w, ignoreNan, dim, 'proportion');

% Calculate the sample proportion (mean)
pHat = compute_mean_helper(x, dim, 'proportion', ignoreNan, w);

% Square the critical z-value for variance calculations
zSquared = zVal .^ 2;

% Calculate the denominator for the Wilson score formula
denominator = 1 + (zSquared ./ n);

% Calculate the adjusted center of the interval
centerAdjusted = pHat + (zSquared ./ (2 .* n));

% Calculate the margin term to add or subtract
margin = zVal .* sqrt((pHat .* (1 - pHat) ./ n) + (zSquared ./ (4 .* (n .^ 2))));

% Apply the margin based on the specified direction and normalize by the denominator
if strcmpi(direction, 'lower')
    bound = (centerAdjusted - margin) ./ denominator;
else
    bound = (centerAdjusted + margin) ./ denominator;
end

% Propagate NaNs if not explicitly ignored
if ~ignoreNan
    bound(any(isnan(x), dim)) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = compute_mean_helper(x, dim, dataType, ignoreNan, w)
% Helper to process conditional mean logic properly

% Compute geometric mean for ratio data by averaging in log domain
if strcmpi(dataType, 'ratio')
    % Log transform for computation
    xLog = log(x);
    
    % Set values that were 0 or negative to NaN since log is undefined there
    xLog(x <= 0) = NaN;
    
    % Compute the standard mean in log space
    yLog = compute_mean_helper(xLog, dim, 'normal', ignoreNan, w);
    
    % Exponentiate back to standard space to get geometric mean
    y = exp(yLog);
    return;
end

% Check if we are dealing with non-normal data or weighted data
if ~strcmpi(dataType, 'normal') || ~isempty(w)
    % Use the extracted stat value from our robust nanstderr function
    y = get_nanstderr_stat_val(x, dim, dataType, w);
elseif ignoreNan
    % Use the built-in nanmean for standard data ignoring NaNs
    y = nanmean(x, dim);
else
    % Use standard mean function
    y = mean(x, dim);
end

% Re-apply NaNs to the result if ignoreNan is false and any NaNs existed
if ~ignoreNan
    y(any(isnan(x), dim)) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = compute_std_helper(x, dim, dataType, ignoreNan, w)
% Helper to route std processing 

% Apply Delta method for standard deviation of ratio data to transform variance cleanly
if strcmpi(dataType, 'ratio')
    % Log transform for computation
    xLog = log(x);
    
    % Set values that were 0 or negative to NaN since log is undefined there
    xLog(x <= 0) = NaN;
    
    % Compute normal mean and normal standard deviation in log space
    meanLog = compute_mean_helper(xLog, dim, 'normal', ignoreNan, w);
    stdLog = compute_std_helper(xLog, dim, 'normal', ignoreNan, w);
    
    % Apply Delta method transformation back to regular space
    y = exp(meanLog) .* stdLog;
    return;
end

% Handle cases where weights are explicitly provided
if ~isempty(w)
    % Create a mask to align valid values between data and weights
    nanMask = isnan(x) | isnan(w);
    xValid = x; 
    wValid = w;
    
    % Apply the NaN mask to both vectors
    xValid(nanMask) = NaN; 
    wValid(nanMask) = NaN;
    
    % Compute the sums of the valid weights and their squares
    wSum = sum(wValid, dim, 'omitnan');
    wSum2 = sum(wValid.^2, dim, 'omitnan');
    
    % Calculate the weighted mean
    mu = sum(wValid .* xValid, dim, 'omitnan') ./ wSum;
    
    % Calculate the unbiased weighted variance
    varW = sum(wValid .* (xValid - mu).^2, dim, 'omitnan') ./ (wSum - (wSum2 ./ wSum));
    
    % Take the square root to get the weighted standard deviation
    y = sqrt(varW);
    
    % Propagate missing values if we are not ignoring them
    if ~ignoreNan
        y(any(isnan(x), dim)) = NaN;
    end
elseif ignoreNan
    % Use the built-in nanstd if we only need to ignore NaNs
    y = nanstd(x, 0, dim);
else
    % Use standard standard deviation function
    y = std(x, 0, dim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = compute_stderr_helper(x, dim, dataType, ignoreNan, w)
% Helper to resolve routing for stderr

% Apply Delta method for standard error of ratio data to transform error cleanly
if strcmpi(dataType, 'ratio')
    % Log transform for computation
    xLog = log(x);
    
    % Set values that were 0 or negative to NaN since log is undefined there
    xLog(x <= 0) = NaN;
    
    % Compute normal mean and normal standard error in log space
    meanLog = compute_mean_helper(xLog, dim, 'normal', ignoreNan, w);
    stderrLog = compute_stderr_helper(xLog, dim, 'normal', ignoreNan, w);
    
    % Apply Delta method transformation back to regular space
    y = exp(meanLog) .* stderrLog;
    return;
end

% Check for complex scenarios requiring nanstderr
if strcmpi(dataType, 'proportion') || ignoreNan || ~isempty(w)
    % Leverage nanstderr for handling weights and proportion specific math
    y = nanstderr(x, dim, 'DataType', dataType, 'Weights', w);
else
    % Use the standard standard error function
    y = stderr(x, dim);
end

% Explicitly propagate NaNs if the user did not opt to ignore them
if ~ignoreNan
    y(any(isnan(x), dim)) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function neff = compute_n_eff(x, w, ignoreNan, dim, dataType)
% Helper to resolve standard or effective sample size computations

% Handle unweighted cases
if isempty(w)
    % Ratios require strictly positive samples
    if strcmpi(dataType, 'ratio')
        neff = sum(x > 0, dim);
    elseif ignoreNan
        % Count elements that are not NaN
        neff = sum(ones(size(x)), dim) - sum(isnan(x), dim);
    else
        % Count all elements along the dimension
        neff = sum(ones(size(x)), dim);
    end
else
    % Identify invalid pairs between data and weights
    nanMask = isnan(x) | isnan(w);
    
    % Mask out invalid weights
    w(nanMask) = NaN;
    
    % Ratios also render non-positive entries invalid
    if strcmpi(dataType, 'ratio')
        w(x <= 0) = NaN;
    end
    
    % Compute sums of the valid weights
    wSum = sum(w, dim, 'omitnan');
    wSum2 = sum(w.^2, dim, 'omitnan');
    
    % Compute Kish's effective sample size
    neff = (wSum.^2) ./ wSum2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = handle_unweighted_only(x, w, statName, unweightedFunc)
% Helper to warn the user and default back when a stat algorithm lacks weight support

% Issue a warning if the user passed weights for an unsupported stat
if ~isempty(w)
    warning('Weights are not currently supported for statistic ''%s''. Ignoring weights.', statName);
end

% Execute the fallback unweighted function
y = unweightedFunc(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = compute_median_helper(x, dim, ignoreNan)
% Wrapper to isolate nanmedian and median behavior for inline handling

% Check if we need to ignore missing values
if ignoreNan
    y = nanmedian(x, dim);
else
    y = median(x, dim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = compute_quartile_helper(x, p, dim, ignoreNan)
% Wrapper to isolate standard quartile behavior vs strict missing-value propagation

% Use standard prctile if ignoring NaNs
if ignoreNan
    y = prctile(x, p, dim);
else
    % Use our strict helper to ensure NaNs propagate
    y = compute_prctile_strict(x, p, dim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = compute_ci_med_helper(x, dim, ignoreNan, statName)
% Helper wrapping the median asymptotic approximation logic

% Compute base stats taking into account whether NaNs are ignored
if ignoreNan
    n = sum(ones(size(x)), dim) - sum(isnan(x), dim);
    med = nanmedian(x, dim);
    iqr = prctile(x, 75, dim) - prctile(x, 25, dim); 
else
    n = sum(ones(size(x)), dim);
    med = median(x, dim);
    iqr = compute_prctile_strict(x, 75, dim) - compute_prctile_strict(x, 25, dim);
end

% Set the constant for the Gaussian-based asymptotic approximation
k = 1.57;

% Apply the approximation based on the specified direction
if strcmpi(statName, 'lower95med')
    y = med - k .* (iqr ./ sqrt(n));
else
    y = med + k .* (iqr ./ sqrt(n));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%