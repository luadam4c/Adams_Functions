function vectorsCell = force_row_cell (vectorsOrig, varargin)
%% Transforms a column cell array or a non-cell array to a row cell array of non-cell vectors
% Usage: vectorsCell = force_row_cell (vectorsOrig, varargin)
% Explanation:
%       This is force_column_cell with 'RowInstead' set to true
%
% Example(s):
%       load_examples;
%       force_row_cell(1:5)
%       force_row_cell(myCellNumeric2D)
%       force_row_cell(myCellColumnVecs)
%       force_row_cell(myCellRowVecs)
%       force_row_cell(myCellStr2D)
%       force_row_cell(myCellStr2D, 'ToLinear', true)
%       force_row_cell(myNumeric2D)
%       force_row_cell(myNumeric3D)
%
% Outputs:
%       vectorsCell - vectors as a row cell array
%                   specified as a row cell array
%
% Arguments:
%       vectorsOrig - original vectors
%                   Note: If a non-cell array, each column is a vector 
%                           to be placed in a cell
%       varargin    - see force_column_cell.m
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%
% Used by:
%       cd/combine_sweeps.m
%       cd/compute_weighted_average.m

% File History:
% 2019-01-13 Modified from force_column_cell.m
% 

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
addRequired(iP, 'vectorsOrig');

% Read from the Input Parser
parse(iP, vectorsOrig, varargin{:});

%% Do the job
vectorsCell = force_column_cell(vectorsOrig, 'RowInstead', true, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
