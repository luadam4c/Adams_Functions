function [valueBoundaries, indBoundaries] = ...
                compute_value_boundaries (values, grouping, varargin)
%% Computes boundaries for values based on a grouping vector
% Usage: [valueBoundaries, indBoundaries] = ...
%               compute_value_boundaries (values, grouping, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       valueBoundaries     - TODO: Description of valueBoundaries
%                           specified as a TODO
%       indBoundaries       - TODO: Description of valueBoundaries
%                           specified as a TODO
%
% Arguments:
%       values      - TODO: Description of reqarg1
%                   must be a TODO
%       grouping    - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for 
%                       compute_index_boundaries()
%
% Requires:
%       cd/compute_index_boundaries.m
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%
% Used by:
%       cd/parse_phase_info.m

% File History:
% 2019-12-02 Created by Adam Lu
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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
% TODO: validation
addRequired(iP, 'values');
addRequired(iP, 'grouping');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, values, grouping, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the compute_index_boundaries() function
otherArguments = iP.Unmatched;

%% Do the job
% Compute index boundaries
indBoundaries = compute_index_boundaries('Grouping', grouping, otherArguments);

% Convert to value units
valueBoundaries = extract_subvectors(values, 'Indices', indBoundaries);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%