function [zScores] = compute_fisher_zscore (corrValues, varargin)
%% Computes the Fisher's z-transformation for an array of correlation values.
% Usage: [zScores] = compute_fisher_zscore(corrValues, varargin)
%
% Explanation:
%       This function transforms correlation coefficients (r) into z-scores (z)
%       using the Fisher z-transformation (arctanh). This transformation
%       stabilizes the variance and makes the sampling distribution of the
%       correlation coefficient approximately normal, which is useful for
%       computing statistics like means and confidence intervals.
%
%       The function can also perform the inverse transformation (tanh) to
%       convert z-scores back into correlation coefficients.
%
% Example(s):
%       % Define some correlation values
%       r = -0.9:0.1:0.9;
%
%       % Perform the forward transformation (r -> z)
%       z = compute_fisher_zscore(r);
%
%       % Perform the inverse transformation (z -> r)
%       r_restored = compute_fisher_zscore(z, 'InverseInstead', true);
%
% Outputs:
%       zScores     - The resulting z-scores (or correlation values if
%                   'InverseInstead' is true). The output will have the
%                   same class (numeric array or cell array) as the input.
%                   specified as a numeric array or a cell array of numeric arrays
%
% Side Effects:
%       None
%
% Arguments:
%       corrValues  - Correlation values (r) or z-scores to be transformed.
%                   For the forward transformation (default), values must
%                   be between -1 and 1, inclusive.
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'InverseInstead': A flag to perform the inverse
%                   transformation (z-scores to correlation values) instead.
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/iscellnumeric.m
%
% Used by:
%       None

% File History:
% 2025-10-15 Created by Gemini
%

%% Hard-coded parameters
% None

%% Default values for optional arguments
inverseInsteadDefault = false;      % default is to perform the forward transform

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
addRequired(iP, 'corrValues', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['First input must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'InverseInstead', inverseInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, corrValues, varargin{:});
inverseInstead = iP.Results.InverseInstead;

%% Check relationships between arguments
% For the forward transformation, correlation values must be within [-1, 1]
if ~inverseInstead
    if iscell(corrValues)
        % Check each cell individually
        cellfun(@(x) validateattributes(x, {'numeric'}, {'>=', -1, '<=', 1}, mfilename, 'corrValues'), corrValues);
    else
        % Check the numeric array
        validateattributes(corrValues, {'numeric'}, {'>=', -1, '<=', 1}, mfilename, 'corrValues');
    end
end

%% Preparation
% Determine if input is a cell array to return output in the same format
isCellInput = iscell(corrValues);

%% Do the job
if inverseInstead
    % Perform the inverse Fisher z-transformation (z -> r)
    if isCellInput
        zScores = cellfun(@tanh, corrValues, 'UniformOutput', false);
    else
        zScores = tanh(corrValues);
    end
else
    % Perform the forward Fisher z-transformation (r -> z)
    if isCellInput
        zScores = cellfun(@atanh, corrValues, 'UniformOutput', false);
    else
        zScores = atanh(corrValues);
    end
end

%% Output results
% The zScores variable is already prepared for output.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%