function objects = plot_ball_stick (radiusSoma, radiusDend, lengthDend, varargin)
%% Plots a ball-and-stick model
% Usage: objects = plot_ball_stick (radiusSoma, radiusDend, lengthDend, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       objects     - rectangle objects plotted
%                   specified as rectangle objects
% Arguments:
%       radiusSoma  - radius of soma
%                   must be a nonegative scalar
%       radiusDend  - radius of dendrite
%                   must be a nonegative scalar
%       lengthDend  - radius of dendrite
%                   must be a nonegative scalar
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2018-11-01 Moved code from TODO
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'radiusSoma', ...
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, radiusSoma, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% Get the current figure
h = gcf;

%% Do the job
% Multiple plots
hold on;

% Plot soma as a circle
objects(1) = rectangle('Position', radiusSoma * [-1, -1, 2, 2], ...
                        'Curvature', [1, 1], 'EdgeColor', edgeColor);

% Plot dendrite as a rectangle
objects(2) = rectangle('Position', [radiusSoma, -radiusDend, ...
                        lengthDend, 2*radiusDend], ...
                        'Curvature', [0, 0], 'EdgeColor', edgeColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%