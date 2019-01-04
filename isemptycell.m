function isEmpty = isemptycell (cellArray, varargin)
%% Returns whether each cell or a cell array is empty
% Usage: isEmpty = isemptycell (cellArray, varargin)
% Explanation:
%       TODO
% Example(s):
%       isemptycell({1:5, [], 2})
%       isemptycell({'', [], 'sdgs'})
% Outputs:
%       isEmpty     - whether each cell or a cell array is empty
%                   specified as a logical array
% Arguments:
%       cellArray   - cell array to test
%                   must be a cell array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/distribute_balls_into_boxes.m
%       cd/find_ind_str_in_cell.m
%       cd/istype.m
%       cd/remove_empty.m

% File History:
% 2010-01-04 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'cellArray', ...
    @(x) validateattributes(x, {'cell'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, cellArray, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
isEmpty = cellfun(@isempty, cellArray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%