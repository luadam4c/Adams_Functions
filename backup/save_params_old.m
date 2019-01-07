function paramTable = save_params (sheetName, paramNames, paramValues, ...
                                paramLowerBounds, paramUpperBounds, varargin)
%% Saves parameters to a spreadsheet file
% Usage: paramTable = save_params (sheetName, paramNames, paramValues, ...
%                               paramLowerBounds, paramUpperBounds, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       paramTable  - a table for all the parameters
%                   specified as a 2-d table
% Arguments:    
%       sheetName   - spreadsheet file name
%                   must be a string scalar or a character vector
%       paramNames  - parameter names
%                   must be a cell vector of character vectors
%       paramValues - parameter values
%                   must be a numeric vector
%       paramLBs    - parameter lower bounds
%                   must be a numeric vector
%       paramUBs    - parameter upper bounds
%                   must be a numeric vector
%
% Requires:
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       TODO cd/isequallength.m
%
% Used by:
%       ~/m3ha/optimizer4gabab/fminsearch3_4compgabab.m
%       ~/m3ha/optimizer4gabab/log_errors_params.m
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

% File History:
% 2018-10-16 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 5
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'sheetName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addRequired(iP, 'paramNames', ...
    @(x) iscellstr(x));
addRequired(iP, 'paramValues', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'paramLowerBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'paramUpperBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, sheetName, paramNames, paramValues, ...
        paramLowerBounds, paramUpperBounds, varargin{:});

% Check relationships between arguments
% TODO: Make a function isequallength.m
if length(unique([numel(paramNames), length(paramValues), ...
                  length(paramLowerBounds), length(paramUpperBounds)])) ~= 1
    fprintf(['Cannot save parameters because names, values, ', ...
             'lower bounds and upper bounds are not all the same length!!\n']);
    paramTable = [];
    return
end

%% Do the job
% Force as a column and rename the variables for the header
Name = force_column_cell(paramNames);
Value = force_column_vector(paramValues);
LowerBound = force_column_vector(paramLowerBounds);
UpperBound = force_column_vector(paramUpperBounds);

% Construct a parameters table
paramTable = table(Name, Value, LowerBound, UpperBound);

% Write the table to the spreadsheet file
writetable(paramTable, sheetName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
