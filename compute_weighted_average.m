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
%       compute_weighted_average([2; NaN; 1], 'Weight', [2, 3, 1], 'AverageMethod', 'rms', 'IgnoreNaN', true)
%       compute_weighted_average(-2, 'AverageMethod', 'rms')
%
% Outputs:
%       avgValues   - averaged value(s)
%                   Note: this is a column vector by default
%                   specified as a numeric vector
%
% Arguments:    
%       values      - values to be averaged (may be a 2D matrix)
%                   must be a numeric array
%       varargin    - 'IgnoreNan': whether to ignore NaN entries
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ForceColumnOutput': whether to force output as a column
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Weights': weights for averaging
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
%       cd/array_fun.m
%       cd/error_unrecognized.m
%       cd/force_row_cell.m
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
% 2019-11-17 Fixed root-mean-square computation
% 2019-11-18 Now returns NaN if values is empty
% 2019-11-18 Fixed 'IgnoreNan' for vectors
% 2019-11-29 Fixed rms average for a single negative value
% 2020-01-02 Added 'ForceColumnOutput' (default == 'true')
% TODO: Simplify math if the weights are all the same 
%       and use this function in compute_stats.m and compute_rms_error.m
% 
tic;

%% Hard-coded parameters
validAverageMethods = {'rms', 'root-mean-square', 'energy', ...
                        'linear', 'arithmetic', 'geometric', 'exponential'};

%% Default values for optional arguments
ignoreNanDefault = false;       % don't ignore NaN by default
forceColumnOutputDefault = true;% force output as a column by default
weightsDefault = [];            % set later
dimToOperateDefault = [];       % use the sum() function default by default
averageMethodDefault = 'rms';   % use the root-mean-square average by default
exponentialWeightingFactorDefault = 50;     % 50%

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
addRequired(iP, 'values', ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IgnoreNan', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceColumnOutput', forceColumnOutputDefault, ...
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
forceColumnOutput = iP.Results.ForceColumnOutput;
valueWeights = iP.Results.Weights;
dimToOperate = iP.Results.DimToOperate;
averageMethod = validatestring(iP.Results.AverageMethod, validAverageMethods);
expWeightingFactor = iP.Results.ExponentialWeightingFactor;

% If averageMethod is 'exponential', any provided weights will be ignored
if strcmpi(averageMethod, 'exponential') && ~isempty(valueWeights)
    message = {['Custom weights are ignored if ', ...
                'averageMethod is ''exponential''!'], ...
                sprintf(['Values will be averaged with ', ...
                        'exponential weighting factor %g%.'], ...
                        expWeightingFactor)};
    mTitle = 'Exponential averaging used';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                            'MessageMode', 'show', 'Verbose', true);
end

%% Preparation
% Return NaN if values is empty
if isempty(values)
    avgValues = NaN;
    return
end

% Remove NaN values if requested
if ignoreNan && ~iscell(values)
    % Test whether each value is NaN
    isNan = isnan(values);

    % Remove NaN values if there is any
    if any(isNan)
        % Otherwise, determine whether to keep each value
        toKeep = isNan;

        % Decide on the dimension to operate
        dimToOperate = decide_on_dimension_to_operate(dimToOperate, values);

        % Align vectors along the dimension to operate 
        if dimToOperate == 1
            % Do nothing
            columnOutput = false;
        elseif dimToOperate == 2
            values = transpose(values);
            toKeep = transpose(toKeep);
            dimToOperate = 1;
            columnOutput = true;
        else
            error('Not implemented yet!');
        end

        if ~any(toKeep)
            % Set the average to be NaN
            if columnOutput || forceColumnOutput
                avgValues = NaN(size(values, 2), 1);
            else
                avgValues = NaN(1, size(values, 2));
            end
            return        
        elseif isvector(values)
            % Values is a vector
            if ~any(toKeep)
                % Set the average to be NaN
                avgValues = NaN;
                return
            else
                % Restrict values to the indices to keep
                values = values(toKeep);

                % Restrict weights to the indices to keep
                if ~isempty(valueWeights)
                    valueWeights = valueWeights(toKeep);
                end
            end
        else
            % Force values as a cell array of column vectors
            values = force_column_vector(values, 'IgnoreNonvectors', false);

            % Force values as a column or row cell array accordingly
            if ~(columnOutput || forceColumnOutput)
                values = force_row_cell(values);
            end

            % Force weights as a column vector
            % TODO: What if weights is a non-vector?
            if ~isempty(valueWeights)
                valueWeights = force_column_vector(valueWeights);
            end
        end
    end
end

%% Do the job
if iscell(values)
    avgValues = ...
        array_fun(@(x) compute_weighted_average(x, ...
                        'IgnoreNan', ignoreNan, ...
                        'Weights', valueWeights, ...
                        'DimToOperate', dimToOperate, ...
                        'AverageMethod', averageMethod, ...
                        'ExponentialWeightingFactor', expWeightingFactor), ...
                values);
else
    avgValues = compute_weighted_average_helper(values, valueWeights, ...
                    dimToOperate, averageMethod, expWeightingFactor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function avgValues = compute_weighted_average_helper (values, valueWeights, ...
                    dimToOperate, averageMethod, expWeightingFactor)

%% Preparation
% If there is only one value, return it
%   Note: This is not true for rms
if numel(values) <= 1
    switch averageMethod
        case {'rms', 'root-mean-square', 'energy'}
            avgValues = abs(values);
        otherwise
            avgValues = values;
    end
    return
end

% Set the dimToOperate that is the default for the sum() function
dimToOperate = decide_on_dimension_to_operate(dimToOperate, values);

% Match up weights to values or simplify things if all weights are the same
if isempty(valueWeights) || numel(unique(valueWeights)) == 1
    valueWeights = [];
    totalWeight = NaN;
else
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
end

%% Do the job
% Compute the average based on various methods
%   See https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/
%           Averaging-types-What-s-the-difference/ta-p/364121
switch averageMethod
    case {'rms', 'root-mean-square', 'energy'}
        % Compute the squared values
        squaredValues = values .* conj(values);

        % Compute the weighted mean of squared values
        if ~isempty(valueWeights)
            meanSquaredValues = sum(valueWeights .* squaredValues, ...
                                        dimToOperate) ./ totalWeight;
        else
            meanSquaredValues = mean(squaredValues, dimToOperate);
        end

        % Compute the weighted root-mean-square average
        avgValues = sqrt(meanSquaredValues);
    case {'linear', 'arithmetic'}
        % Compute the weighted linear average
        if ~isempty(valueWeights)
            avgValues = sum(valueWeights .* values, dimToOperate) ./ ...
                            totalWeight;
        else
            avgValues = mean(values, dimToOperate);
        end
    case {'geometric'}
        % Compute the weighted geometric average
        if ~isempty(valueWeights)
            avgValues = exp(sum(valueWeights .* log(values), dimToOperate) ...
                            ./ totalWeight);
        else
            avgValues = exp(mean(log(values), dimToOperate));
        end
    case {'exponential'}
        % Compute the reciprocal of temperature T 
        %   Note: T == 1/(1-EWF), so w = 1/T == 1-EWF
        w = 1 - (expWeightingFactor / 100);

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

function dimToOperate = decide_on_dimension_to_operate (dimToOperate, values)
%% Set the dimToOperate that is the default for the sum() function

if isempty(dimToOperate)
    % This is the first dimension that is not one
    dimToOperate = find(size(values) ~= 1, 1, 'first');
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

% Set default weights
if isempty(valueWeights)
    valueWeights = ones(size(values));
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
