function vectorsCell = match_format_cell (vectorsOrig, cellToMatch, varargin)
%% Matches the format of a set of vectors to a cell array to match so that cellfun can be used
% Usage: vectorsCell = match_format_cell (vectorsOrig, cellToMatch, varargin)
% Explanation:
%       This function is a combination of force_column_cell.m
%           and match_dimensions.m
% Example(s):
%       TODO
% Outputs:
%       vectorsCell - vectors as a cell array with matching dimensions
%                   specified as a column cell array
% Arguments:
%       vectorsOrig - original vectors
%                   Note: If an array, each column is a vector 
%                           to be placed in a cell
%                   must be a numeric array or a cell array 
%                       or a character vector
%       vectorsOrig - a cell array whose dimensions are to be matched
%                   must be a cell array
%
% Requires:
%       cd/force_column_cell.m
%       cd/match_dimensions.m
%
% Used by:    

% File History:
% 2018-10-31 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

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

% Add required inputs to the Input Parser
addRequired(iP, 'vectorsOrig', ...
    @(x) isnumeric(x) || iscell(x) || ischar(x));
addRequired(iP, 'cellToMatch', ...
    @(x) iscell(x));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, vectorsOrig, cellToMatch, varargin{:});

%% Do the job
% Make sure vectorsOrig is a column cell array
vectorsCell = force_column_cell(vectorsOrig);

% Match the dimensions with cellToMatch
vectorsCell = match_dimensions(vectorsCell, size(cellToMatch));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Make sure numeric vectors are column vectors
% Note: this is now done in force_column_cell.m as well
if isnumeric(vectorsOrig)
    force_column_numeric(vectorsOrig);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%