function avgValues = compute_weighted_average (values, varargin)
%% Computes a weighted average value (root-mean-square by default)
% Usage: avgValues = compute_weighted_average (values, varargin)
% Explanation:
%   See https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/
%           Averaging-types-What-s-the-difference/ta-p/364121
%           for the explanation on different averaging methods
% Example(s):
%       TODO
% Outputs:
%       avgValues   - averaged value(s)
%                   specified as a numeric vector
%
% Arguments:    
%       values      - values to be averaged (may be a 2D matrix)
%                   must be a numeric array
%       varargin    - 'Weights': weights
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
%                       'exponential' - exponential averaging
%                   default == 'rms'
%                   - 'ExponentialWeightingFactor': weighting factor (%)
%                           for the exponential method
%                   must be a numeric scalar between 0~100
%                   default == 50 %
%
% Requires:
%       cd/error_unrecognized.m
%       cd/get_variable_name.m
%       cd/ispositiveintegerscalar.m
%       cd/match_dimensions.m
%
% Used by:    
%       cd/compute_sweep_errors.m
%
% Related functions:
%       cd/compute_means.m
%       cd/compute_rms_error.m

% File History:
% 2018-10-26 Created by Adam Lu
% 2018-10-28 Fixed the case when values has less than one element
% TODO: Simply math if the weights are all the same 
%       and use this function in compute_means.m and compute_rms_error.m
% 

%% Hard-coded parameters
validAverageMethods = {'rms', 'root-mean-square', 'energy', ...
                        'linear', 'arithmetic', 'exponential'};

%% Default values for optional arguments
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
    @(x) validateattributes(x, {'numeric'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
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
    case {'exponential'}
        % Compute the reciprocal of temperature T 
        %   Note: T == 1/(1-EWF), so w = 1/T == 1-EWF
        w = 1 - exponentialWeightingFactor;

        % Compute the weighted exponential average
        % TODO: deal with dimension not 1 and ndims not 2
        % TODO: make into a function
        for iRow = 1:nRows
            if iRow == 1
                avgValues = values(iRow, :);
            else
                avgValues = w * values(iRow, :) + (1 - w) * avgValues;
            end
        end
    otherwise
        error_unrecognized(get_variable_name(averageMethod), ...
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