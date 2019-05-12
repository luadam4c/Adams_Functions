function stats = compute_stats (vecs, statName, varargin)
%% Computes a statistic of vector(s) possibly restricted by endpoint(s)
% Usage: stats = compute_stats (vecs, statName, dim (opt), varargin)
% Explanation:
%       TODO
% Example(s):
%       data = randn(10, 3);
%       compute_stats(data, 'mean')
%       compute_stats(data, 'std')
%       compute_stats(data, 'stderr')
%       compute_stats(data, 'lower95')
%       compute_stats(data, 'upper95')
%       compute_stats(data, 'cov')
%       compute_stats(data, 'zscore')
%       compute_stats(data, 'mean', 2)
% Outputs:
%       stats       - the computed statistic for each vector
%                   specified as a numeric vector 
% Arguments:
%       vecs        - vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%       statName    - name of the statistic
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'mean'      - mean
%                       'std'       - standard deviation
%                       'stderr'    - standard error
%                       'lower95'   - lower bound of the 95% confidence interval
%                       'upper95'   - upper bound of the 95% confidence interval
%                       'cov'       - coefficient of variation
%                       'zscore'    - z-score
%                       'range'     - range
%                       'range2mean'- range over mean
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
%       cd/remove_outliers.m
%       cd/stderr.m
%       cd/nanstderr.m
%
% Used by:
%       cd/parse_multiunit.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
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
% TODO: Combine with compute_weighted_average.m
% 

%% Hard-coded parameters
validStatNames = {'mean', 'std', 'stderr', 'lower95', 'upper95', ...
                    'cov', 'zscore', 'range', 'range2mean'};

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
    case 'mean'
        if ignoreNan
            func = @(x) nanmean(x, dim);
        else
            func = @(x) mean(x, dim);
        end            
    case 'std'
        if ignoreNan
            func = @(x) nanstd(x, dim);
        else
            func = @(x) std(x, dim);
        end            
    case 'stderr'
        if ignoreNan
            func = @(x) nanstderr(x, dim);
        else
            func = @(x) stderr(x, dim);
        end            
    case 'lower95'
        if ignoreNan
            func = @(x) nanmean(x, dim) - 1.95 .* nanstderr(x, dim);
        else
            func = @(x) mean(x, dim) - 1.95 .* stderr(x, dim);
        end            
    case 'upper95'
        if ignoreNan
            func = @(x) nanmean(x, dim) + 1.95 .* nanstderr(x, dim);
        else
            func = @(x) mean(x, dim) + 1.95 .* stderr(x, dim);
        end     
    case 'cov'
        if ignoreNan
            func = @(x) nanstd(x, dim) ./ nanmean(x, dim);
        else
            func = @(x) std(x, dim) ./ mean(x, dim);
        end
    case 'zscore'
        if ignoreNan
            func = @(x) nanmean(x, dim) ./ nanstd(x, dim);
        else
            func = @(x) mean(x, dim) ./ std(x, dim);
        end
    case 'range'
        func = @(x) range(x, dim);
    case 'range2mean'
        if ignoreNan
            func = @(x) range(x, dim) ./ nanmean(x, dim);
        else
            func = @(x) range(x, dim) ./ mean(x, dim);
        end
    otherwise
        error('Code logic error!');
end

% Compute the statistic for each subvector
if iscell(subVecs)
    stats = cellfun(func, subVecs);
else
    stats = func(subVecs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
