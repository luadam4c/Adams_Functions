function fig = decide_on_fighandle (varargin)
%% Decides on the figure handle depending on what's provided
% Usage: fig = decide_on_fighandle (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       fig         - figure handle to use
%                   specified as a Figure object handle
% Arguments:
%       varargin    - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - Any other parameter-value pair for the figure() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_bar.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-05-10 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
figHandleDefault = [];              % no existing figure by default
figNumberDefault = [];              % no figure number by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));

% Read from the Input Parser
parse(iP, varargin{:});
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;

% Keep unmatched arguments for the figure() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
if ~isempty(figHandle)
    % Use the given figure
    fig = figure(figHandle, otherArguments{:});
elseif ~isempty(figNumber)
    % Create and clear a new figure with given figure number
    fig = figure(figNumber, otherArguments{:});
    clf(fig);
else
    % Get the current figure or create one if non-existent
    fig = gcf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%