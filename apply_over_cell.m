function output = apply_over_cell (myFunction, cellArray, varargin)
%% Apply a function that usually takes two equivalent arguments over all contents of a cell array
% Usage: output = apply_over_cell (myFunction, cellArray, varargin)
% Explanation:
%       TODO
% Examples:
%       vecs1 = {[2, 3, 4], [3, 4, 5], [1, 3, 4]};
%       apply_over_cell(@intersect, vecs1)
%       apply_over_cell(@union, vecs1, 'OptArg', 'stable')
%       load_examples;
%       apply_over_cell(@outerjoin, myCellTable, 'MergeKeys', true);
% Outputs:
%       output      - final output
%
% Arguments:
%       myFunction  - a custom function that takes two equivalent arguments
%                       as normal input
%                       e.g., intersect(), union(), outerjoin()
%                   must be a function handle
%       cellArray   - a cell array of arguments for myFunction
%                   must be a cell array of inputs that 
%                       can serve as the first two inputs for myFunction
%       varargin    - 'OptArg': optional argument that's 
%                                   not a parameter-value pair 
%                   default == []
%                   - optional parameter-value pairs for myFunction
%
%
% Used by:
%       cd/combine_variables_across_tables.m
%       ~/m3ha/optimizer4gabab/compare_and_plot_across_conditions.m
%
% File History:
% 2019-03-17 MOdified from union_over_cell.m

%% Hard-coded parameters

%% Default values for optional arguments
optArgDefault = '';         % Not provided

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'myFunction', ...           % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
addRequired(iP, 'cellArray', @iscell);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OptArg', optArgDefault);

% Read from the Input Parser
parse(iP, myFunction, cellArray, varargin{:});
optArg = iP.Results.OptArg;

% Keep unmatched arguments for myFunction
otherArguments = struct2arglist(iP.Unmatched);

%% Return the input if not a cell array
if isempty(cellArray) || ~iscell(cellArray)
    output = cellArray;
    return;
end

%% Return the union of the contents of the cell array
% Count the number of elements
nElements = numel(cellArray);

% If there are no elements, return an empty matrix
if nElements == 0
    output = [];
    return
end

% Initialize the output
output = cellArray{1};

% Iterate over all elements
if nElements > 1
    for iElement = 2:nElements
        if ~isempty(optArg)
            output = myFunction(output, cellArray{iElement}, ...
                                optArg, otherArguments{:});
        else
            output = myFunction(output, cellArray{iElement}, otherArguments{:});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%