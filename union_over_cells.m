function unionOfContents = union_over_cells (cellArray, varargin)
%% Apply the union function over all contents of a cell array
% Usage: unionOfContents = union_over_cells (cellArray, varargin)
% Explanation:
%       TODO
%
% Examples:
%       union_over_cells([2; 4; 4])
%       union_over_cells({'ab', 'bc', 'bg'})
%       union_over_cells({[1 2], [2; 3]})
%
% Outputs:    
%       unionOfContents - union of contents
%
% Arguments:
%       cellArray   - a cell array of arrays that will be unioned;
%                       if just an array, return the array
%                   must be a cell array of input arrays that 
%                       can be recognized by the built-in union function
%       varargin    - 'SetOrder': 
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'sorted' - sort the values
%                       'stable' - use the original order of the values
%                   default == 'sorted'
%                   - Any other parameter-value pair for the union() function
%
%
% Used by:
%       cd/create_default_grouping.m
%       cd/plot_tuning_curve.m
%       ~/m3ha/optimizer4gabab/compare_and_plot_across_conditions.m

% File History:
%   2018-08-17 Created by Adam Lu
%   2019-08-21 Now uses apply_over_cells.m

%% Hard-coded parameters
validSetOrders = {'sorted', 'stable'};

%% Default values for optional arguments
setOrderDefault = 'sorted';

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'cellArray');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SetOrder', setOrderDefault, ...
    @(x) any(validatestring(x, validSetOrders)));

% Read from the Input Parser
parse(iP, cellArray, varargin{:});
setOrder = validatestring(iP.Results.SetOrder, validSetOrders);

% Keep unmatched arguments for the union() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
unionOfContents = apply_over_cells(@union, cellArray, 'OptArg', setOrder, ...
                                    otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
