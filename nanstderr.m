function varargout = nanstderr(X, varargin)
%% Calculate the standard error of the mean or proportion excluding NaN values
% Usage: [stdErr, statVal] = nanstderr(X, varargin)
% Explanation:
%       Calculates the standard error of the mean or proportion for elements 
%       in a numeric array X, treating NaNs as missing values. Optionally 
%       incorporates a weight vector for weighted standard error calculations.
%
% Example(s):
%       stdErr = nanstderr([1, 2, NaN, 4])
%       [stdErr, statVal] = nanstderr([1, 0, NaN, 1], 'DataType', 'proportion')
%       [stdErr, statVal] = nanstderr([1, 0, 1, 1], 'DataType', 'proportion', 'Weights', [0.5, 1, 1, 2])
%
% Outputs:
%       stdErr      - standard error excluding NaN values
%                   specified as a numeric array
%       statVal     - (opt) calculated mean or proportion (weighted if specified)
%                   specified as a numeric array
%
% Side Effects:
%       None
%
% Arguments:
%       X           - input data
%                   must be a numeric array
%       varargin    - 'dim': dimension to operate along
%                   must be a positive numeric vector
%                   default == first non-singleton dimension of X
%                   - 'DataType': type of data to calculate standard error for
%                   must be an unambiguous, case-insensitive match to one of:
%                       'normal'        - standard error of the mean
%                       'proportion'    - standard error of a proportion
%                   default == 'normal'
%                   - 'Weights': optional weight vector for calculating 
%                       weighted standard errors
%                   must be empty or a numeric array of the same size as X
%                   default == []
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_stats.m
%       cd/plot_grouped_jitter.m
%       cd/ZG_fit_IEI_distributions.m
%       /media/adamX/Paula_IEIs/paula_iei4.m
%       scAAV/analyze_qupath_figure2.m
%       scAAV/analyze_qupath_figure3.m

% File History:
% 2017-12-14 Created
% 2019-03-14 Added 'dim' as an optional argument
% 2025-09-11 Fixed the case where dim is 2
% 2026-03-03 Fixed the case where dim is 2
% 2026-03-03 Added 'DataType' and 'Weights' optional arguments
% 2026-03-03 Added statVal as a second output argument
% 2026-03-04 Now uses mean(), std(), and sum() with the 'omitnan' option

%% Hard-coded parameters
validDataTypes = {'normal', 'proportion'};

%% Default values for optional arguments
dimDefault = [];
dataTypeDefault = 'normal';
weightsDefault = [];

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
addRequired(iP, 'X', ...                   % vectors
    @(x) assert(isnum(x), 'X must be a numeric array!'));

% Add optional inputs to the Input Parser
addOptional(iP, 'dim', dimDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataType', dataTypeDefault, ...
    @(x) any(validatestring(x, validDataTypes)));
addParameter(iP, 'Weights', weightsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x), 'Weights must be empty or a numeric array!'));

% Read from the Input Parser
parse(iP, X, varargin{:});
dim = iP.Results.dim;
dataType = validatestring(iP.Results.DataType, validDataTypes);
weights = iP.Results.Weights;

%% Do the job
if isempty(dim)
    % Find the first non-singleton dimension, matching standard nanstd behavior
    dim = find(size(X) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

if isempty(weights)
    % Calculate number of valid observations
    nValid = sum(~isnan(X), dim);

    % Calculate standard error and conditionally the statistic based on data type
    switch dataType
        case 'normal'
            % Calculate standard error of the mean
            stdErr = std(X, 0, dim, 'omitnan') ./ sqrt(nValid);
            
            % Compute the mean only if requested
            if nargout > 1
                statVal = mean(X, dim, 'omitnan');
            end
        case 'proportion'
            % Calculate the proportion
            p = mean(X, dim, 'omitnan');
            
            % Calculate standard error of the proportion
            stdErr = sqrt(p .* (1 - p) ./ nValid);
            
            % Store the proportion if requested
            if nargout > 1
                statVal = p;
            end
    end
else
    % Validate weights are same size as X or allow implicit expansion
    assert(isequal(size(weights), size(X)) || isvector(weights), ...
        'Weights must be the same size as X or a vector!');
        
    % Mask NaNs consistently across X and weights
    nanMask = isnan(X) | isnan(weights);
    X(nanMask) = NaN;
    weights(nanMask) = NaN;
    
    % Calculate the sum of weights and sum of squared weights
    wSum = sum(weights, dim, 'omitnan');
    wSum2 = sum(weights.^2, dim, 'omitnan');
    
    % Calculate standard error and conditionally the statistic based on data type
    switch dataType
        case 'normal'
            % Calculate the weighted mean
            mu = sum(weights .* X, dim, 'omitnan') ./ wSum;

            % Calculate the weighted sample variance
            varW = sum(weights .* (X - mu).^2, dim, 'omitnan') ./ (wSum - (wSum2 ./ wSum));

            % Calculate Kish's effective sample size
            nEff = (wSum.^2) ./ wSum2;

            % Calculate the standard error of the weighted mean
            stdErr = sqrt(varW ./ nEff);
            
            % Store the weighted mean if requested
            if nargout > 1
                statVal = mu;
            end
        case 'proportion'
            % Calculate the weighted proportion
            p = sum(weights .* X, dim, 'omitnan') ./ wSum;

            % Calculate the standard error of the weighted proportion
            stdErr = sqrt(p .* (1 - p) .* wSum2 ./ (wSum.^2));
            
            % Store the weighted proportion if requested
            if nargout > 1
                statVal = p;
            end
    end
end

% Assign outputs
varargout{1} = stdErr;
if nargout > 1
    varargout{2} = statVal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%