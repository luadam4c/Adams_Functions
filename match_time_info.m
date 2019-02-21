function [tVecs, siMs] = match_time_info (tVecs, siMs, varargin)
%% Match time vector(s) and sampling interval(s)
% Usage: [tVecs, siMs] = match_time_info (tVecs, siMs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       tVecs       - time vector(s)
%                   specified as a numeric array 
%                       or a cell array of numeric arrays
%       siMs        - sampling interval(s) in milliseconds
%                   specified as a positive vector
% Arguments:
%       tVecs       - time vector(s)
%                   must be empty or 
%                       a numeric array or a cell array of numeric arrays
%       siMs        - sampling interval(s) in milliseconds
%                   must be empty or a positive vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/compute_sampling_interval.m
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/match_row_count.m
%       cd/match_vector_count.m
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-02-20 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'tVecs');
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either empty, a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs');
    @(x) assert(isempty(x) || ispositivevector(x), ...
                ['siMs must be either empty or a positive vector!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, tVecs, siMs, varargin{:});
param1 = iP.Results.param1;

%% Do the job
% Compute sampling interval(s) and create time vector(s)
if isempty(siMs) && ~isempty(tVecs)
    % Compute sampling interval(s) in ms
    siMs = compute_sampling_interval(tVecs);
elseif isempty(tVecs) && ~isempty(siMs)
    % Create time vector(s)
    tVecs = create_time_vectors(nSamples, 'SamplingIntervalMs', siMs, ...
                                    'TimeUnits', 'ms');
elseif isempty(tVecs) && isempty(siMs)
    error('One of siMs and tVecs must be provided!');
end

% Count the maximum number of entries
nEntries = max(count_vectors(tVecs), numel(siMs));

% Match the number of entries
siMs = match_row_count(siMs, nEntries);
tVecs = match_vector_count(tVecs, nEntries);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%