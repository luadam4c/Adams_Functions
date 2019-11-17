function avgValues = compute_weighted_average (values, varargin)
%% Computes a weighted average value (root-mean-square by default)
% Usage: avgValues = compute_weighted_average (values, varargin)
% Explanation:
%   See https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/
%           Averaging-types-What-s-the-difference/ta-p/364121
%           for the explanation on different averaging methods
%
% Example(s):
%       compute_weighted_average([1; 10; 100], 'AverageMethod', 'rms')
%       compute_weighted_average([1; 10; 100], 'AverageMethod', 'linear')
%       compute_weighted_average([1; 10; 100], 'AverageMethod', 'geometric')
%       compute_weighted_average([1; 10; 100], 'AverageMethod', 'exponential')
%       compute_weighted_average([100; 10; 1], 'AverageMethod', 'exponential')
%       compute_weighted_average([NaN, 3, 27; NaN, 4, 64], 'AverageMethod', 'geometric', 'DimToOperate', 2, 'IgnoreNan', true)
%       compute_weighted_average([NaN; 10], 'Weight', [2, 1], 'AverageMethod', 'linear', 'IgnoreNaN', true)
%
% Outputs:
%       avgValues   - averaged value(s)
%                   specified as a numeric vector
%
% Arguments:    
%       values      - values to be averaged (may be a 2D matrix)
%                   must be a numeric array
%       varargin    - 'IgnoreNan': whether to ignore NaN entries
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Weights': weights
%                   must be a numeric array
%                   default == ones(size(values))
%                   - 'DimToOperate': dimension to compute averages on
%                   must be empty or a positive integer scalar
%                   default == [] (uses the sum() function default)
%                   - 'AverageMethod': method for averaging
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'rms', 'root-mean-square', 'energy'
%                                     - root-mean-square averaging
%                       'linear', 'arithmetic' 
%                                     - linear averaging
%                       'geometric'   - geometric averaging
%                       'exponential' - exponential averaging
%                   default == 'rms'
%                   - 'ExponentialWeightingFactor': weighting factor (%)
%                           for the exponential method
%                   must be a numeric scalar between 0~100
%                   default == 50 %
%
% Requires:
%       cd/error_unrecognized.m
%       cd/force_column_vector.m
%       cd/get_var_name.m
%       cd/ispositiveintegerscalar.m
%       cd/match_dimensions.m
%
% Used by:
%       cd/test_normality.m
%       cd/compute_lts_errors.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%
% Related functions:
%       cd/compute_stats.m
%       cd/compute_rms_error.m

% File History:
% 2018-10-26 Created by Adam Lu
% 2018-10-28 Fixed the case when values has less than one element
% 2019-01-11 Added 'geometric' as an averaging method
% 2019-10-10 Added 'IgnoreNan' as an optional argument
% 2019-10-12 Allow values to be a cell array
% 2019-10-12 Fixed 'IgnoreNan' for matrices
% TODO: Simply math if the weights are all the same 
%       and use this function in compute_stats.m and compute_rms_error.m
% 

%% Hard-coded parameters
validAverageMethods = {'rms', 'root-mean-square', 'energy', ...
                        'linear', 'arithmetic', 'geometric', 'exponential'};

%% Default values for optional arguments
ignoreNanDefault = false;       % don't ignore NaN by default
weightsDefault = [];            % set later
dimToOperateDefault = [];       % use the sum() function default by default
averageMethodDefault = 'rms';   % use the root-mean-square average by default
exponentialWeightingFactorDefault = 50;     % 50%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'values', ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IgnoreNan', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Weights', weightsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'3d'}));
addParameter(iP, 'DimToOperate', dimToOperateDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['DimToOperate must be either empty ', ...
                    'or a positive integer scalar!']));

addParameter(iP, 'AverageMethod', averageMethodDefault, ...
    @(x) any(validatestring(x, validAverageMethods)));
addParameter(iP, 'ExponentialWeightingFactor', ...
                exponentialWeightingFactorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 100}));

% Read from the Input Parser
parse(iP, values, varargin{:});
ignoreNan = iP.Results.IgnoreNan;
valueWeights = iP.Results.Weights;
dimToOperate = iP.Results.DimToOperate;
averageMethod = validatestring(iP.Results.AverageMethod, validAverageMethods);
exponentialWeightingFactor = iP.Results.ExponentialWeightingFactor;

% If averageMethod is 'exponential', any provided weights will be ignored
if strcmpi(averageMethod, 'exponential') && ~isempty(valueWeights)
    message = {['Custom weights are ignored if ', ...
                'averageMethod is ''exponential''!'], ...
                sprintf(['Values will be averaged with ', ...
                        'exponential weighting factor %g%.'], ...
                        exponentialWeightingFactor)};
    mTitle = 'Exponential averaging used';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                            'MessageMode', 'show', 'Verbose', true);
end

%% Preparation
% Remove NaN values if requested
if ignoreNan && ~iscell(values)
    if isvector(values)
        if all(isnan(values))
            avgValues = NaN;
            return
        else
            values = values(~isnan(values));
            if ~isempty(valueWeights)
                valueWeights = valueWeights(~isnan(values));
            end
        end
    else
        % Force as a cell array of column vectors
        if dimToOperate == 1
            values = force_column_vector(values, 'IgnoreNonvectors', false);
        elseif dimToOperate == 2
            values = force_column_vector(transpose(values), ...
                                            'IgnoreNonvectors', false);
        else
            error('Not implemented yet!');
        end
        if ~isempty(valueWeights)
            valueWeights = force_column_vector(valueWeights);
        end
        dimToOperate = 1;
    end
end

%% Do the job
if iscell(values)
    avgValues = ...
        cellfun(@(x) compute_weighted_average(x, ...
                'IgnoreNan', ignoreNan, ...
                'Weights', valueWeights, ...
                'DimToOperate', dimToOperate, ...
                'AverageMethod', averageMethod, ...
                'ExponentialWeightingFactor', exponentialWeightingFactor), ...
            values);
else
    avgValues = compute_weighted_average_helper(values, valueWeights, ...
                    dimToOperate, averageMethod, exponentialWeightingFactor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function avgValues = compute_weighted_average_helper (values, valueWeights, ...
                    dimToOperate, averageMethod, exponentialWeightingFactor)

%% Preparation
% If there is only one value, return it
if numel(values) <= 1
    avgValues = values;
    return
end

% Set default weights
if isempty(valueWeights)
    valueWeights = ones(size(values));
end

% Set the dimToOperate that is the default for the sum() function
if isempty(dimToOperate)
    % This is the first dimension that is not one
    dimToOperate = find(size(values) ~= 1, 1, 'first');
end

% If valueWeights is a vector, make sure it is along the dimension to operate
if isvector(valueWeights)
    if dimToOperate == 1
        valueWeights = reshape(valueWeights, [length(valueWeights), 1, 1]);
    elseif dimToOperate == 2
        valueWeights = reshape(valueWeights, [1, length(valueWeights), 1]);
    elseif dimToOperate == 3
        valueWeights = reshape(valueWeights, [1, 1, length(valueWeights)]);
    end
end

% Make sure the dimensions of the weights array 
%   does not differ from the values array in the dimension to operate
if size(valueWeights, dimToOperate) ~= size(values, dimToOperate)
    error(['The dimensions of the weights array cannot differ ', ...
            'from that of the values array in the dimension to operate!']);    
end

% Match up the weights array to the dimensions of the values array
valueWeights = match_dimensions(valueWeights, size(values));

% Sum all the weights
totalWeight = sum(valueWeights, dimToOperate);

%% Do the job
% Compute the average based on various methods
%   See https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/
%           Averaging-types-What-s-the-difference/ta-p/364121
switch averageMethod
    case {'rms', 'root-mean-square', 'energy'}
        % Compute the weighted root-mean-square average
        avgValues = sqrt(sum(valueWeights .* values .^ 2, dimToOperate) ./ ...
                        totalWeight);
    case {'linear', 'arithmetic'}
        % Compute the weighted linear average
        avgValues = sum(valueWeights .* values, dimToOperate) ./ totalWeight;
    case {'geometric'}
        % Compute the weighted geometric average
        avgValues = exp(sum(valueWeights .* log(values), dimToOperate) ...
                        ./ totalWeight);
    case {'exponential'}
        % Compute the reciprocal of temperature T 
        %   Note: T == 1/(1-EWF), so w = 1/T == 1-EWF
        w = 1 - (exponentialWeightingFactor / 100);

        % Compute the weighted exponential average
        % TODO: deal with dimension not 1 and ndims not 2
        % TODO: make into a function
        for iRow = 1:size(values, 1)
            if iRow == 1
                avgValues = values(iRow, :);
            else
                avgValues = w * values(iRow, :) + (1 - w) * avgValues;
            end
        end
    otherwise
        error_unrecognized(get_var_name(averageMethod), ...
                            averageMethod, mfilename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Make sure the weights array and the values array 
%   have an equal number of dimensions
if ndims(valueWeights) ~= ndims(values)
    error(['The weights array must have an equal number of ', ...
            'dimensions as the values array!']);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
