function [outTables] = combine_variables_across_tables (inTables, varargin)
%% Combines measures across different tables
% Usage: [outTables] = combine_variables_across_tables (inTables, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       outTables   - tables organized with each 
%                   specified as a TODO
% Arguments:
%       inTables    - tables organized with each measure as a column
%                   must be a cell array of tables
%       varargin    - 'VariableNames': variable (column) names of the table
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == plot all variables
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/clc2_plot_measures.m

% File History:
% 2019-03-15 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
variableNamesDefault = {};  % plot all variables by default

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
addRequired(iP, 'inTables', ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VariableNames', variableNamesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));

% Read from the Input Parser
parse(iP, inTables, varargin{:});
varNames = iP.Results.VariableNames;

%% Preparation
% TODO

%% Do the job
% TODO

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%