function [output1] = plot_horizontal_line (yValue, varargin)
%% Plots a horizontal line
% Usage: [output1] = plot_horizontal_line (yValue, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       yValue     - TODO: Description of yValue
%                   must be a TODO
%       varargin    - 'XLimits': TODO: Description of XLimits
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 201X-XX-XX Created by TODO or Adapted from TODO
% TODO: Use unMatched varargin parts as parameters for line()
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of XLimits

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
addRequired(iP, 'yValue', ...                  % TODO: Description of yValue
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XLimits', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, yValue, varargin{:});
xLimits = iP.Results.XLimits;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

% Check relationships between arguments
% TODO

%% Preparation
% TODO

%% Do the job
% TODO

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%