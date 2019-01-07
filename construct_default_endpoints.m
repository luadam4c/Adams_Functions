function endPoints = construct_default_endpoints (nSamples, varargin)
%% Constructs default endpoints from number of samples
% Usage: endPoints = construct_default_endpoints (nSamples, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       endPoints   - end points for each vector
%                   specified as a numeric array 
%                       or a cell array of numeric arrays
% Arguments:
%       nSamples    - number of samples for each vector
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/iscellnumeric.m
%
% Used by:
%       cd/extract_subvectors.m

% File History:
% 2019-01-03 Created by Adam Lu
% TODO: Make 'StartIndex' and optional argument with default == 1
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
addRequired(iP, 'nSamples', ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['nSamples must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, nSamples, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if isempty(nSamples)
    endPoints = [];
elseif iscell(nSamples)
    endPoints = cellfun(@(x) construct_default_endpoints(x), ...
                        nSamples, 'UniformOutput', false);
else
    endPoints = transpose([ones(size(nSamples)), nSamples]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%