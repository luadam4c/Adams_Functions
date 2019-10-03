function grouping = create_grouping_by_vectors (array, varargin)
%% Creates a grouping array that matches input by putting all elements of a vector into the same group
% Usage: grouping = create_grouping_by_vectors (array, varargin)
% Explanation:
%       Creates an array with the same dimensions as the input array
%           but with each entry replaced by the vector number
%
% Example(s):
%       grouping1 = create_grouping_by_vectors(magic(5))
%       grouping2 = create_grouping_by_vectors({1:5, 2:3, 6:10})
%       grouping3 = create_grouping_by_vectors({1:5, 1:2; 1:3, 1:4})
%       grouping4 = create_grouping_by_vectors({{1:5}, {1:3, 1:4}})
%
% Outputs:
%       grouping    - a grouping array with the same size as the input array
%                   specified as a numeric array
%
% Arguments:
%       array       - an array
%       varargin    - 'TreatCellAsArray': whether to treat a cell array
%                                           as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellNumAsArray': whether to treat a cell array
%                                       of numeric arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellStrAsArray': whether to treat a cell array
%                                       of character arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/iscellnumeric.m
%       cd/isnum.m
%
% Used by:
%       cd/create_default_grouping.m
%       cd/plot_swd_histogram.m

% File History:
% 2018-12-27 Created by Adam Lu
% 2019-01-15 Improved code performance
% 2019-10-03 Now allows array to be cell arrays
% 

%% Hard-coded parameters

%% Default values for optional arguments
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellNumAsArrayDefault = false; % treat cell arrays of numeric arrays
                                    %   as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default

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
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellNumAsArray', treatCellNumAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, array, varargin{:});
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;

%% Do the job
if isnum(array) || iscell(array) && treatCellAsArray || ...
        iscellnumeric(array) && treatCellNumAsArray || ...
        iscellstr(array) && treatCellStrAsArray
    % Count the number of rows and columns
    nRows = size(array, 1);
    nColumns = size(array, 2);

    % Use the column number
    grouping = ones(nRows, 1) * (1:nColumns);
elseif iscellnumeric(array)
    % Count the number of vectors
    nVectors = numel(array);

    % Create an array of group numbers
    groupNumbers = reshape(1:nVectors, size(array));


    % Give each vector a different group number
    grouping = cellfun(@(x, y) y * ones(size(x)), ...
                        array, num2cell(groupNumbers), 'UniformOutput', false);
elseif iscell(array)
    % Do this for each cell
    grouping = cellfun(@(x) create_grouping_by_vectors(x, ...
                            'TreatCellAsArray', treatCellAsArray, ...
                            'TreatCellNumAsArray', treatCellNumAsArray, ...
                            'TreatCellStrAsArray', treatCellStrAsArray), ...
                        array, 'UniformOutput', false);
else
    error('Input has no vectors!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% This is two times slower:
grouping = repmat(1:nColumns, [nRows, 1])

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%