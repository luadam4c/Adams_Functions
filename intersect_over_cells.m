function intersection = intersect_over_cells (cellArray, varargin)
%% Apply the intersect function over all contents of a cell array
% Usage: intersection = intersect_over_cells (cellArray, varargin)
% Explanation:
%       TODO
%
% Examples:
%       intersect_over_cells([2; 3])
%       intersect_over_cells({'ab', 'bc', 'bg'})
%       intersect_over_cells({[1 2], [2; 3]})
%
% Outputs:    
%       intersection    - intersection of contents
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
%                   - Any other parameter-value pair for intersect()
%
% Used by:
%       TODO

% File History:
%   2018-01-10 Created
%   2018-05-08 Fixed the case if cellArray is not a cell array
%   2018-05-08 Changed tabs to spaces and limited width to 80
%   2018-08-17 Renamed intersectm -> intersect_over_cells
%   2019-08-21 Now uses apply_over_cells.m

%% Hard-coded parameters
validSetOrders = {'sorted', 'stable'};

%% Default values for optional arguments
setOrderDefault = 'sorted';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
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

% Keep unmatched arguments for the intersect() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
intersection = apply_over_cells(@intersect, cellArray, 'OptArg', setOrder, ...
                                    otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
