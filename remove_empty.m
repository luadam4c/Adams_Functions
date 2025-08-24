function array = remove_empty (array, varargin)
%% Removes empty elements from an array
% Usage: array = remove_empty (array, varargin)
% Explanation:
%       Removes elements considered 'empty' from an array
%           For numeric arrays NaN is considered 'empty'
%
% Example(s):
%       remove_empty([NaN 1 2 NaN])
%       remove_empty({'a', '', 'b'})
%       remove_empty({{'a', '', 'b'}, {}, {'c', 'd'}})
%
% Outputs:
%       array       - cleaned array
% Arguments:
%       array       - original array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isemptycell.m
%
% Used by:
%       cd/extract_elements.m
%       cd/parse_all_abfs.m
%       cd/unique_custom.m

% File History:
% 2019-01-04 Created by Adam Lu
% 2025-08-23 Fixed the condition that array is not a cell array
% TODO: Allow option to replace by NaNs
% TODO: Somehow preserve matrix structure if an array is not a vector
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

%% Do the job
if iscell(array)
    array = array(~isemptycell(array));
else
    array = array(~isnan(array));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%