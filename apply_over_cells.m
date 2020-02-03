function output = apply_over_cells (myFunction, inputs, varargin)
%% Apply a function that usually takes two equivalent arguments over all contents of a cell array
% Usage: output = apply_over_cells (myFunction, inputs, varargin)
% Explanation:
%       TODO
%
% Examples:
%       vecs1 = {[2, 3, 4], [3, 4, 5], [1, 3, 4]};
%       apply_over_cells(@intersect, vecs1)
%       apply_over_cells(@union, vecs1, 'OptArg', 'stable')
%       load_examples;
%       apply_over_cells(@outerjoin, myCellTable, 'MergeKeys', true)
%       apply_over_cells(@outerjoin, myCellTable, 'MergeKeys', true, 'Keys', 'Key')
%
% Outputs:
%       output      - final output
%
% Arguments:
%       myFunction  - a custom function that takes two equivalent arguments
%                       as normal input
%                       e.g., intersect(), union(), outerjoin()
%                   must be a function handle
%       inputs      - a cell array of arguments for myFunction
%                   must be a cell array of inputs that 
%                       can serve as the first two inputs for myFunction
%       varargin    - 'OptArg': optional argument that's 
%                                   not a parameter-value pair 
%                   default == []
%                   - optional parameter-value pairs for myFunction
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/combine_variables_across_tables.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_rank_neurons.m
%       cd/plot_vertical_line.m
%       ~/m3ha/optimizer4gabab/m3ha_compare_and_plot_across_conditions.m.m
%
% File History:
%   2019-03-17 Modified from union_over_cell.m
%   2019-08-21 Now accepts any array and 
%               applies it to itself if not a cell array

%% Hard-coded parameters

%% Default values for optional arguments
optArgDefault = '';         % Not provided

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'myFunction', ...           % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
addRequired(iP, 'inputs');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OptArg', optArgDefault);

% Read from the Input Parser
parse(iP, myFunction, inputs, varargin{:});
optArg = iP.Results.OptArg;

% Keep unmatched arguments for myFunction
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Return the input if empty
if isempty(inputs)
    output = inputs;
    return;
end

% Count the number of elements
nElements = numel(inputs);

% If there are no elements, return an empty matrix
if nElements == 0
    output = [];
    return
end

% If inputs is not a cell array, apply the function to itself
if ~iscell(inputs)
    if ~isempty(optArg)
        output = myFunction(inputs, inputs, optArg, otherArguments{:});
    else
        output = myFunction(inputs, inputs, otherArguments{:});
    end
    return
end

%% Iterate over all cells
% Initialize the output as the first input
output = inputs{1};

% Iterate over all elements
if nElements > 1
    for iElement = 2:nElements
        if ~isempty(optArg)
            output = myFunction(output, inputs{iElement}, ...
                                optArg, otherArguments{:});
        else
            output = myFunction(output, inputs{iElement}, otherArguments{:});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
