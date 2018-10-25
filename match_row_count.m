function arrayNew = match_row_count (arrayOld, nRowsNew, varargin)
%% Expands or truncates an array to match a given number of rows (dimension #1)
% Usage: arrayNew = match_row_count (arrayOld, nRowsNew, varargin)
% Explanation:
%       TODO
% Example(s):
%       match_row_count([1, 2, 3], 6)
%       match_row_count([1; 2; 3], 6)
%       match_row_count([1, 2; 3, 4], 7)
% Outputs:
%       arrayNew    - array matched
%                   specified as a numeric, cell or struct array
%
% Arguments:    
%       arrayOld    - array to match
%                   must be a numeric, cell or struct array
%       nRowsNew    - new number of rows
%                   must be a positive integer scalar
%
% Used by:    
%       cd/m3ha_compute_single_neuron_errors.m

% File History:
% 2018-10-25 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'arrayOld', ...
    @(x) validateattributes(x, {'numeric', 'cell', 'struct'}, {'3d'}));
addRequired(iP, 'nRowsNew', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, arrayOld, nRowsNew, varargin{:});

%% Preparation
% Query the old number of rows
nRowsOld = size(arrayOld, 1);

% If the new number of rows are the same as the old ones, 
%   just return the old array
if nRowsNew == nRowsOld
    arrayNew = arrayOld;
    return
end

% Query the number of dimensions
nDims = ndims(arrayOld);

%% Expand or truncate array
if nRowsNew > nRowsOld
    % Compute the factor to expand by
    factorToExpand = floor(nRowsNew / nRowsOld);

    % Compute the remaining number of rows to fill after expansion
    remainingNrows = mod(nRowsNew, nRowsOld);

    % Expand array
    if nDims == 2
        % First expand by factorToExpand
        arrayNew = repmat(arrayOld, [factorToExpand, 1]);

        % Fill in remaining rows by the first rows
        arrayNew = vertcat(arrayNew, arrayOld(1:remainingNrows, :));
    elseif nDims == 3
        % First expand by factorToExpand
        arrayNew = repmat(arrayOld, [factorToExpand, 1, 1]);

        % Fill in remaining rows by the first rows
        arrayNew = vertcat(arrayNew, arrayOld(1:remainingNrows, :, :));
    end
elseif nRowsNew < nRowsOld
    % Truncate array
    if nDims == 2
        arrayNew = arrayOld(nRowsNew, :);
    elseif nDims == 3
        arrayNew = arrayOld(nRowsNew, :, :);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Query the dimensions
dims = size(arrayOld);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%