function [combArray, paramsUsed] = ...
                compute_combined_array (arrays, combineMethod, varargin)
%% Computes a combined array from a cell array of numeric arrays
% Usage: [combArray, paramsUsed] = ...
%               compute_combined_array (arrays, combineMethod, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       arrs = {magic(3), magic(3), transpose(magic(3))};
%       compute_combined_array(arrs, 'mean')
%       compute_combined_array(arrs, 'average')
%       compute_combined_array(arrs, 'std')
%       compute_combined_array(arrs, 'stderr')
%       compute_combined_array(arrs, 'lower95')
%       compute_combined_array(arrs, 'upper95')
%       compute_combined_array(arrs, 'maximum')
%       compute_combined_array(arrs, 'minimum')
%       compute_combined_array(arrs, 'sum')
%       compute_combined_array(arrs, 'prod')
%       compute_combined_array(arrs, 'all')
%       compute_combined_array(arrs, 'any')
%
% Outputs:
%       combArray   - the combined array
%                   specified as a numeric array
%
% Arguments:
%       arrays      - arrays to average
%                   must be a numeric array or a cell array or numeric arrays
%       combineMethod   - method for combining arrays
%                       must be an unambiguous, case-insensitive match to one of: 
%                           'average' or 'mean' - take the average
%                           'std'       - standard deviation
%                           'stderr'    - standard error
%                           'lower95'   - lower bound of the 95% conf interval
%                           'upper95'   - upper bound of the 95% conf interval
%                           'maximum'   - take the maximum
%                           'minimum'   - take the minimum
%                           'sum'       - take the sum
%                           'prod'      - take the product
%                           'all'       - take the logical AND
%                           'any'       - take the logical OR
%       varargin    - 'IgnoreNaN': whether to ignore NaN values
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/apply_over_cells.m
%       cd/compute_stats.m
%       cd/create_error_for_nargin.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2021-05-23 Created by Adam Lu
% 

%% Hard-coded parameters
validCombineMethods = {'average', 'mean', ...
                        'std', 'stderr', 'lower95', 'upper95', ...
                        'maximum', 'minimum', 'sum', 'prod', ...
                        'all', 'any'};

%% Default values for optional arguments
ignoreNanDefault = true;            % ignore NaN values by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'arrays', ...
    @(x) assert(isnum(x) || iscell(x), ...
                'arrays must be either a numeric array or a cell array!'));
addRequired(iP, 'combineMethod', ...
    @(x) any(validatestring(x, validCombineMethods)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IgnoreNaN', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));

% Read from the Input Parser
parse(iP, arrays, combineMethod, varargin{:});
ignoreNan = iP.Results.IgnoreNaN;

% Validate combine method
combineMethod = validatestring(combineMethod, validCombineMethods);

%% Preparation
% Make sure arrays is a cell array
if ~iscell(arrays)
    arrays = {arrays};
end

% Find the number of dimensions of the array
nDimsAll = cellfun(@ndims, arrays);
maxNDims = max(nDimsAll);

%% Do the job
% Decide on the dimension number to concatenate arrays
dimToCat = maxNDims + 1;

% Concatenate all arrays
catted = apply_over_cells(@(x, y) cat(dimToCat, x, y), arrays);

% Combine arrays
switch combineMethod
    case {'average', 'mean'}
        % Take the mean of each row and return a column
        combArray = compute_stats(catted, 'mean', dimToCat, ...
                                    'IgnoreNan', ignoreNan);
    case {'std', 'stderr', 'lower95', 'upper95'}
        combArray = compute_stats(catted, combineMethod, dimToCat, ...
                                    'IgnoreNan', ignoreNan);
    case {'maximum', 'minimum'}
        if ignoreNan
            nanFlag = 'omitnan';
        else
            nanFlag = 'includenan';
        end
        
        switch combineMethod
            case 'maximum'
                % Take the maximum of all elements across catted dimension
                combArray = max(catted, [], dimToCat, nanFlag);
            case 'minimum'
                % Take the minimum of all elements across catted dimension
                combArray = min(catted, [], dimToCat, nanFlag);
        end
    case 'sum'
        % Add all elements across catted dimension
        combArray = sum(catted, dimToCat);
    case 'prod'
        % Multiply all elements across catted dimension
        combArray = prod(catted, dimToCat);
    case 'all'
        % Take the logical AND of all elements across catted dimension
        combArray = all(catted, dimToCat);
    case 'any'
        % Take the logical OR of all elements across catted dimension
        combArray = any(catted, dimToCat);
end

%% Output info
paramsUsed.combineMethod = combineMethod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%