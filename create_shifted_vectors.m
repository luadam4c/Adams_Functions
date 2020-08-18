function vecsNew = create_shifted_vectors (vecOrig, shiftValues, varargin)
%% Creates shifted vectors based on values to shift
% Usage: vecsNew = create_shifted_vectors (vecOrig, shiftValues, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       create_shifted_vectors(-2:2, 3)
%
% Outputs:
%       vecsNew     - TODO: Description of vecsNew
%                   specified as a TODO
%
% Arguments:
%       vecOrig     - TODO: Description of vecOrig
%                   must be a TODO
%       shiftValues - TODO: Description of vecOrig
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%
% Used by:
%       cd/compute_relative_event_times.m

% File History:
% 2020-08-18 Moved from compute_relative_event_times.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'vecOrig');
addRequired(iP, 'shiftValues');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vecOrig, shiftValues, varargin{:});
% param1 = iP.Results.param1;


%% Do the job
% Force the original vector as a column vector
vecOrig = force_column_vector(vecOrig);

% Force the values to shift as a row vector
shiftValues = force_row_vector(shiftValues);

% Create shifted vectors
%   Note: Each column corresponds to a time window
vecsNew = repmat(shiftValues, size(vecOrig)) + ...
            repmat(vecOrig, size(shiftValues));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%