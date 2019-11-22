function [avgTrace, paramsUsed] = compute_maximum_trace (traces, varargin)
%% Computes the maximum of traces that are not necessarily the same length
% Usage: [avgTrace, paramsUsed] = compute_maximum_trace (traces, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       avgTrace    - the maximum trace
% Arguments:    
%       traces      - traces to maximum
%                   Note: If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%       varargin    - Any other parameter-value pair for compute_combined_trace()
%                   
% Requires:
%       cd/compute_combined_trace.m
%
% Used by:
%       cd/create_indices.m

% File History:
% 2019-01-03 Adapted from compute_average_trace.m
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
addRequired(iP, 'traces');

% Read from the Input Parser
parse(iP, traces, varargin{:});

% Keep unmatched arguments for compute_combined_trace()
otherArguments = iP.Unmatched;

%% Do the job
[avgTrace, paramsUsed] = ...
    compute_combined_trace(traces, 'maximum', otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
