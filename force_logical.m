function arrayNew = force_logical (arrayOld, varargin)
%% Forces any numeric binary array to become a logical array
% Usage: arrayNew = force_logical (arrayOld, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       arrayNew    - new array
%                   specified as an array
% Arguments:
%       arrayOld    - old array
%                   must be an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/isbinaryarray.m
%
% Used by:
%       cd/load_params.m

% File History:
% 2018-12-11 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'arrayOld');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, arrayOld, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if isnumeric(arrayOld) && isbinaryarray(arrayOld)
    % Force as logical
    arrayNew = logical(arrayOld);
else
    % Return original array
    arrayNew = arrayOld;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%