function set_visible_off (object)
%% Set the 'Visible' property of object(s) to 'off'
% Usage: set_visible_off (object)
% Explanation:
%       TODO
%
% Example(s):
%       yAxises = get(ax, 'YAxis');
%       set_visible_off(yAxises);
%
% Arguments:
%       object  - a MATLAB object
%               must be an object with a 'Visible' property
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_traces.m
%       cd/update_figure_for_corel.m

% File History:
% 2019-12-22 Moved from update_figure_for_corel.m
% 

%% Hard-coded parameters

%% Default values for optional arguments

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
addRequired(iP, 'object');

% Read from the Input Parser
parse(iP, object);

%% Do the job
if iscell(object)
    cellfun(@set_visible_off, object);
else
    set(object, 'Visible', 'off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%