function column = force_column (array, varargin)
%% Forces an array as a column
% Usage: column = force_column (array, varargin)
% Explanation:
%       Forces an array as a column so that match_row_count.m can be used
%       Calls force_column_cell.m or force_column_numeric.m 
%           depending on the type of the array
% Example(s):
%       TODO
% Outputs:
%       column      - array reorganized as a column
%                   specified as a column vector
% Arguments:
%       array       - array to reorganize
%                   must be an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
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
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'array', ...                  % TODO: Description of array
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, array, varargin{:});
param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

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