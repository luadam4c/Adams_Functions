function figHandle = update_figure_for_corel (figHandle, varargin)
%% Update figure to be journal-friendly
% Usage: figHandle = update_figure_for_corel (figHandle, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       fig = update_figure_for_corel(fig);
%
% Outputs:
%       figHandle     - TODO: Description of figHandle
%                   specified as a TODO
%
% Arguments:
%       figHandle     - TODO: Description of figHandle
%                   must be a TODO
%       varargin    - 'RemoveTicks': whether to remove all ticks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/struct2arglist.m
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-19 Created by Adam Lu
% 2019-09-20 Added 'RemoveTicks' as an optional argument
% 

%% Hard-coded parameters
% TODO: Make optional parameters
fontSizeLabels = 8;
fontSizeAxis = 6; 
linewidth = 1;
tickLengths = [0.025, 0.025];

%% Default values for optional arguments
removeTicksDefault = false;  % set later

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
addRequired(iP, 'figHandle');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RemoveTicks', removeTicksDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, figHandle, varargin{:});
removeTicks = iP.Results.RemoveTicks;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
% TODO

%% Preparation
% TODO

titleFontSizeMultiplier = fontSizeLabels / fontSizeAxis;
labelFontSizeMultiplier = fontSizeLabels / fontSizeAxis;

%% Do the job
% Find all axes in the figure
ax = findall(figHandle, 'type', 'axes');

% Count the number of axes
nAx = numel(ax);

% Set font
set(ax, 'FontName', 'Arial');
set(ax, 'FontSize', fontSizeAxis);
set(ax, 'TitleFontSizeMultiplier', 4/3);
set(ax, 'TitleFontWeight', 'normal');
set(ax, 'LabelFontSizeMultiplier', 4/3);

% Make all rulers linewidth 1
set(ax, 'LineWidth', 1);

% Remove boxes
for iAx = 1:nAx
    box(ax(iAx), 'off');
end

% Make all ticks go outward
set(ax, 'TickDir', 'out');
set(ax, 'TickDirMode', 'manual');
set(ax, 'TickLength', tickLengths);

% Remove x and y axis ticks
if removeTicks
    set(ax, 'XTick', []);
    set(ax, 'YTick', []);
end

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

for iAx = 1:nAx
    ax(iAx).XAxis.LineWidth = 1;
    ax(iAx).YAxis.LineWidth = 1;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%