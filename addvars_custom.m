function myTable = addvars_custom (myTable, value, varargin)
%% Adds a column to a table, matching rows if necessary 
% Usage: myTable = addvars_custom (myTable, value, varargin)
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
%       cd/m3ha_plot_figure07.m

% File History:
% 2020-03-06 Created by Adam Lu

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

% Add the column to the table
myTable = addvars(myTable, value, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%