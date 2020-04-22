function stats = compute_stats (vecs, statName, varargin)
%% Computes a statistic of vector(s) possibly restricted by endpoint(s)
% Usage: stats = compute_stats (vecs, statName, dim (opt), varargin)
% Explanation:
%       TODO
%       Note: If any element is empty, returns NaN
%
% Example(s):
%       data = randn(10, 3);
%       compute_stats(data, 'mean')
%       compute_stats(data, 'std')
%       compute_stats(data, 'stderr')
%       compute_stats(data, 'err')
%       compute_stats(data, 'lower95')
%       compute_stats(data, 'upper95')
%       compute_stats(data, 'cov')
%       compute_stats(data, 'zscore')
%       compute_stats(data, 'mean', 2)
%       compute_stats(data, 'max', 2)
%       compute_stats(data, 'mean', 'IgnoreNan', true)
%
% Outputs:
%       stats       - the computed statistic for each vector
%                   specified as a numeric vector 
% Arguments:
%       vecs        - vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%       statName    - name of the statistic
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'average' or 'mean' - mean
%                       'std'       - standard deviation
%                       'stderr'    - standard error
%                       'err'       - error margin
%                       'lower95'   - lower bound of the 95% confidence interval
%                       'upper95'   - upper bound of the 95% confidence interval
%                       'cov'       - coefficient of variation
%                       'zscore'    - z-score
%                       'max'       - maximum
%                       'min'       - manimum
%                       'range'     - range
%                       'range2mean'- percentage of range relative to mean
%       dim         - (opt) dimension to compute stats along
%                   must be either 1, 2 or 3
%                   default == 1
%       varargin    - 'IgnoreNan': whether to ignore NaN entries
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
%       cd/select_similar_values.m
%       cd/test_difference.m
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
% TODO: Combine with compute_weighted_average.m
% 

%% Hard-coded parameters
validStatNames = {'average', 'mean', 'std', 'stderr', 'err', ...
                    'lower95', 'upper95', ...
                    'cov', 'zscore', 'max', 'min', 'range', 'range2mean'};

%% Default values for optional arguments
dimDefault = 1;                 % compute across rows by default
ignoreNanDefault = false;       % don't ignore NaN by default
removeOutliersDefault = false;  % don't remove outliers by default
indicesDefault = [];            % set later
endPointsDefault = [];          % set later
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
addRequired(iP, 'statName', ...
    @(x) any(validatestring(x, validStatNames)));

% Add optional inputs to the Input Parser
addOptional(iP, 'dim', dimDefault, ...
    @(x) assert(isnumeric(x) && (x == 1 || x == 2 || x == 3), ...
                'dim must be either 1, 2 or 3!'));

% Add parameter-value pairs to the Input Parser
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
dim = iP.Results.dim;
ignoreNan = iP.Results.IgnoreNan;
removeOutliers = iP.Results.RemoveOutliers;
indices = iP.Results.Indices;
endPoints = iP.Results.EndPoints;
windows = iP.Results.Windows;

%% Do the job
% Extract subvectors
subVecs = extract_subvectors(vecs, 'Indices', indices, ...
                            'EndPoints', endPoints, 'Windows', windows);

% Remove outliers if requested
if removeOutliers
    if iscell(subVecs)
        subVecs = cellfun(@remove_outliers, subVecs, ...
                        'UniformOutput', false);
    else
        subVecs = remove_outliers(subVecs);
    end
end

% Decide on the function to use on each vector
switch statName
    case {'average', 'mean'}
        if ignoreNan
            func = @(x) nanmean(x, dim);
        else
            func = @(x) mean(x, dim);
        end            
    case 'std'
        if ignoreNan
            func = @(x) nanstd(x, 0, dim);
        else
            func = @(x) std(x, 0, dim);
        end            
    case 'stderr'
        if ignoreNan
            func = @(x) nanstderr(x, dim);
        else
            func = @(x) stderr(x, dim);
        end
    case {'err', 'lower95', 'upper95'}
        % Compute the number of samples along the dimension
        if ignoreNan
            nFunc = @(x) sum(ones(size(x)), dim) - sum(isnan(x), dim);
        else
            nFunc = @(x) sum(ones(size(x)), dim);
        end

        % Compute the multiple of standard errors for the t distribution 95% 
        %   confidence interval
        tFunc = @(x) arrayfun(@(y) tinv(0.975, y), nFunc(x));

        switch statName
        case 'err'
            if ignoreNan
                func = @(x) tFunc(x) .* nanstderr(x, dim);
            else
                func = @(x) tFunc(x) .* stderr(x, dim);
            end
        case 'lower95'
            if ignoreNan
                func = @(x) nanmean(x, dim) - tFunc(x) .* nanstderr(x, dim);
            else
                func = @(x) mean(x, dim) - tFunc(x) .* stderr(x, dim);
            end
        case 'upper95'
            if ignoreNan
                func = @(x) nanmean(x, dim) + tFunc(x) .* nanstderr(x, dim);
            else
                func = @(x) mean(x, dim) + tFunc(x) .* stderr(x, dim);
            end
        end
    case 'cov'
        if ignoreNan
            func = @(x) nanstd(x, 0, dim) ./ nanmean(x, dim);
        else
            func = @(x) std(x, 0, dim) ./ mean(x, dim);
        end
    case 'zscore'
        if ignoreNan
            func = @(x) nanmean(x, dim) ./ nanstd(x, 0, dim);
        else
            func = @(x) mean(x, dim) ./ std(x, 0, dim);
        end
    case 'max'
        func = @(x) max(x, [], dim);
    case 'min'
        func = @(x) min(x, [], dim);
    case 'range'
        func = @(x) range(x, dim);
    case 'range2mean'
        if ignoreNan
            func = @(x) (range(x, dim) ./ nanmean(x, dim)) * 100;
        else
            func = @(x) (range(x, dim) ./ mean(x, dim)) * 100;
        end
    otherwise
        error('Code logic error!');
end

% Compute the statistic for each subvector
if iscell(subVecs)
    % If any element is empty, return with NaN
    if ~any(isemptycell(subVecs))
        stats = cellfun(func, subVecs);
    else
        stats = nan(size(subVecs));
        parfor iVec = 1:numel(subVecs)
            if ~isempty(subVecs{iVec})
                stats(iVec) = func(subVecs{iVec});
            end
        end
    end
else
    if isempty(subVecs)
        stats = NaN;
    else
        stats = func(subVecs);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

case 'lower95'
    if ignoreNan
        func = @(x) nanmean(x, dim) - 1.96 .* nanstderr(x, dim);
    else
        func = @(x) mean(x, dim) - 1.96 .* stderr(x, dim);
    end            
case 'upper95'
    if ignoreNan
        func = @(x) nanmean(x, dim) + 1.96 .* nanstderr(x, dim);
    else
        func = @(x) mean(x, dim) + 1.96 .* stderr(x, dim);
    end     

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
