function myTable = addvar_as_rowname (myTable, value, varargin)
%% Adds a column to a table in the beginning and as row name
% Usage: myTable = addvar_as_rowname (myTable, value, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       myTable     - a table with column added
%                   specified as a table
%
% Arguments:
%       myTable     - a table
%                   must be a table
%       value       - value(s) for the column to add
%                   must be a scalar or a vector
%       varargin    - Any other parameter-value pair for addvars()
%
% Requires:
%       cd/match_row_count.m
%
% Used by:

% File History:
% 2020-08-18 Adapted from addvars_custom.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

%% Do the job
% Count the number of rows
nRows = height(myTable);

% Make sure value has the same number of rows
value = match_row_count(value, nRows);

% Add as row names
myTable.Properties.RowNames = value;

% Add column to the table
myTable = addvars(myTable, value, 'NewVariableNames', inputname(2), ...
                    'Before', 1, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
