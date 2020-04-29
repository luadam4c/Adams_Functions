function elements1 = match_positions (array1, array2, elements2, varargin)
%% Finds element(s) of an array that matches the positions of elements in a second list
% Usage: elements1 = match_positions (array1, array2, elements2, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       match_positions({45, 15, 2}, {'cars', 'dogs', 'bat'}, {'o', 'a', 's'})
%       label = match_positions(labels, types, type);
%       tau1 = match_positions(coeffValues, coeffNames, 'b');
%       tau2 = match_positions(coeffValues, coeffNames, 'd');
%
% Outputs:
%       elements1   - matched elements
%
% Arguments:
%       array1      - an array
%                   must be an array
%       array2      - a second array
%                   must be an array
%       elements2   - elements in the second array
%                   must be an array
%       varargin    - 'MaxNum': maximum number of positions to match
%                   must be a positive integer scalar or Inf
%                   default == Inf
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/find_first_match.m
%       cd/ispositiveintegerscalar.m
%
% Used by:
%       cd/compute_peak_decay.m
%       cd/extract_channel.m
%       cd/m3ha_network_launch.m
%       cd/plot_table_parallel.m
%       cd/update_figure_for_corel.m

% File History:
% 2018-12-15 Created by Adam Lu
% 2018-12-24 Now accepts any array type as the first element
% 2019-12-30 Now accepts any array type as the second and third arguments
% TODO: Use interp1_custom if only numeric arrays passed in
%       Note: Change the argument order to be consistent with interp1
%               Should be array2, array1, elements2

%% Default values for optional arguments
maxNumDefault = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'array1');
addRequired(iP, 'array2');
addRequired(iP, 'elements2');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MaxNum', maxNumDefault, ...
    @(x) isinf(x) || ispositiveintegerscalar(x));

% Read from the Input Parser
parse(iP, array1, array2, elements2);
maxNum = iP.Results.MaxNum;

%% Do the job
% For each element in elements2, find the index of the first match in array2
idxMatch = find_first_match(elements2, array2);

% Get all the elements of array1 in this position
elements1 = extract_subvectors(array1, 'Indices', idxMatch);

% Restrict to the maximum number of elements
if numel(array1) > maxNum
    elements1 = array1(1:maxNum);
end

% If a cell array, return as an element if there is only one element
if iscell(elements1) && numel(elements1) == 1
    elements1 = elements1{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%