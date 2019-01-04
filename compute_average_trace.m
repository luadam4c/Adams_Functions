function [avgTrace, paramsUsed] = compute_average_trace (traces, varargin)
%% Computes the average of traces that are not necessarily the same length
% Usage: [avgTrace, paramsUsed] = compute_average_trace (traces, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       avgTrace    - the average trace
% Arguments:    
%       traces      - traces to average
%                   Note: If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%       varargin    - Any other parameter-value pair for compute_combined_trace()
%                   
% Requires:
%       cd/compute_combined_trace.m
%
% Used by:
%       cd/compute_average_pulse_response.m
%       cd/find_passive_params.m
%       cd/force_column_numeric.m
%       cd/m3ha_import_raw_traces.m

% File History:
% 2018-10-11 Created by Adam Lu
% 2018-12-18 Now uses extract_subvectors.m
% 2019-01-03 Moved code to compute_combined_trace.m
% 

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'traces');

% Read from the Input Parser
parse(iP, traces, varargin{:});

% Keep unmatched arguments for compute_combined_trace()
otherArguments = iP.Unmatched;

%% Do the job
[avgTrace, paramsUsed] = ...
    compute_combined_trace(traces, 'average', otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
