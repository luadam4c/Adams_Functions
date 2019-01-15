function grouping = create_grouping_by_columns (array, varargin)
%% Creates a grouping array using column numbers
% Usage: grouping = create_grouping_by_columns (array, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       grouping    - a grouping array with the same size as the input array
%                   specified as a numeric array
% Arguments:
%       array       - an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/create_default_grouping.m
%       cd/plot_swd_histogram.m

% File History:
% 2018-12-27 Created by Adam Lu
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
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, array, varargin{:});
% param1 = iP.Results.param1;

%% Preparation
% Count the number of rows and columns
nRows = size(array, 1);
nColumns = size(array, 2);

%% Do the job
% Use the column number
grouping = repmat(1:nColumns, [nRows, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%