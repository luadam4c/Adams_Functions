function nAsEachC = count_A_each_C (nAsEachB, nBsEachC)
%% Counts the number of (A)s in each (C) based on the number of (A)s in each (B) and the number of (B)s in each (C)
% Usage: nAsEachC = count_A_each_C (nAsEachB, nBsEachC)
% Explanation:
%       TODO
% Example(s):
%       count_A_each_C([3, 4, 3, 2, 1], [2; 3])
%       count_A_each_C([3, 4, 3, 2, 1, 4], [2; 3])
% Outputs:
%       nAsEachC    - number of (A)s in each (C)
%                   specified as a numeric vector
% Arguments:
%       nAsEachB    - number of (A)s in each (B)
%                   must be  a numeric vector
%       nBsEachC    - the number of (B)s in each (C)
%                   must be  a numeric vector
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/combine_data_from_same_slice.m

% File History:
% 2019-08-21 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'nAsEachB', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'nBsEachC', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, nAsEachB, nBsEachC);

%% Preparation
% Check if the number of (B)s match up
if length(nAsEachB) ~= sum(nBsEachC)
    fprintf('The number of (B)s provided do not match up!\n');
    nAsEachC = [];
    return
end

%% Do the job
% Get the indices of all the last (B)s for each (C)
iLastBEachC = cumsum(nBsEachC);

% Get the indices of all the first (B)s for each (C)
iFirstBEachC = iLastBEachC - nBsEachC + 1;

% Count the number of (A)s in each (C)
nAsEachC = arrayfun(@(x, y) sum(nAsEachB(x:y)), iFirstBEachC, iLastBEachC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%