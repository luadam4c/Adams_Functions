function [isInTable] = is_var_in_table (varName, table, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: [isInTable] = is_var_in_table (varName, table, varargin)
% Explanation:
%       TODO
% Example(s):
%       T = readtable('outages.csv');
%       is_var_in_table('Cause', T)
%
% Outputs:
%       isInTable     - TODO: Description of isInTable
%                   specified as a TODO
% Arguments:
%       varName     - TODO: Description of varName
%                   must be a TODO
%       table       - TODO: Description of table
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/struct2arglist.m
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 201X-XX-XX Created by TODO or Adapted from TODO
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'varName');
addRequired(iP, 'table');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, varName, table, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
% TODO

%% Preparation
% TODO

%% Do the job
% TODO

allVarNames = T.Properties.VariableName
str= allVarNames
ismatch(varName,allVarNames)
A = [0.53 0.67 0.01 0.38 0.07 0.42 0.69];
B = any()
 
% data type: character arrays
%A character array is a sequence of characters
% C= 'hello, world'
%convert different type of data to character array: 
%C= char (A)

% data type: cell
%cell array is a data type with indexed data containers called cells, 
%where each cell can contain any type of data. 
%Cell arrays commonly contain either lists of text, 
%combinations of text and numbers, or numeric arrays of different sizes.
%C= {x, y, z}

% data type: table
%table arrays store column-oriented or tabular data, 
%such as columns from a text file or spreadsheet. 
%Tables store each piece of column-oriented data in a variable.

% T.Properties.VariableNames= {variable1, v2, v3,...}

% ismatch()
%Returns whether each element in a list matches a candidate
%Usage: [isMatch, indices, matched] = ismatch(list, cand, varargin)

% any()
%True if any element of a vector is a nonzero number or is
%    logical 1 (TRUE).  any ignores entries that are NaN (Not a Number).
%  For vectors, any(V) returns logical 1 (TRUE) if any of the 
%  elements of the vector is a nonzero number or is logical 1 (TRUE).
%  Otherwise it returns logical 0 (FALSE).
%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%