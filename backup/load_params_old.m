function [paramNames, paramValues, paramLBs, paramUBs] = ...
                read_params (sheetName, varargin)
%% Loads parameters to a spreadsheet file
% Usage: [paramNames, paramValues, paramLBs, paramUBs] = ...
%               read_params (sheetName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       paramNames  - parameter names
%                   specified as a column cell vector of character vectors
%       paramValues - parameter values
%                   specified as a column numeric vector
%       paramLBs    - parameter lower bounds
%                   specified as a column numeric vector
%       paramUBs    - parameter upper bounds
%                   specified as a column numeric vector
% Arguments:    
%       sheetName   - spreadsheet file name
%                   must be a string scalar or a character vector
%
% Used by:    
%       /TODO:dir/TODO:file

% File History:
% 2018-10-16 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

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

% Add required inputs to the Input Parser
addRequired(iP, 'sheetName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, sheetName, varargin{:});

%% Do the job
% Read the spreadsheet file into a table
neuronParamsTable = readtable(sheetName);

% Read from the table
paramNames = neuronParamsTable.Name;
paramValues = neuronParamsTable.Value;
paramLBs = neuronParamsTable.LowerBound;
paramUBs = neuronParamsTable.UpperBound;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
