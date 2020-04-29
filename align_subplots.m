function figHandle = align_subplots (figHandle, varargin)
%% Aligns subplots in a figure
% Usage: figHandle = align_subplots (figHandle (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       figHandle   - handle to updated figure
%                   specified as a Figure object handle
%
% Arguments:
%       figHandle   - (opt) figure handle to update
%                   must be a Figure object handle
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for create_subplots()
%
% Requires:
%       cd/argfun.m
%       cd/create_subplots.m
%       cd/extract_elements.m
%       cd/remove_non_axes.m
%
% Used by:
%       cd/update_figure_for_corel.m

% File History:
% 2019-12-29 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
figHandleDefault = [];
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'figHandle', figHandleDefault);

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, figHandle, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the create_subplots() function
otherArguments = iP.Unmatched;

%% Do the job
% Find all axes in the figure
ax = findall(figHandle, 'type', 'axes');

% Remove titles and labels
ax = remove_non_axes(ax);

% Count the number of axes
nAx = numel(ax);

% Get all axes outer positions (this returns a cell array)
axOuterPositions = get(ax, 'OuterPosition');

% Extract the lefts, bottoms, widths and heights
[outerLefts, outerBottoms, outerWidths, outerHeights] = ...
    argfun(@(x) extract_elements(axOuterPositions, 'specific', 'Index', x), ...
            1, 2, 3, 4);

% Combine into a table
posTable = table(ax, outerLefts, outerBottoms, outerWidths, outerHeights);

% Sort by the bottoms in descending order, then by the lefts in ascending order
posTable = sortrows(posTable, [3, 2], {'descend', 'ascend'});

% Extract fields
ax = posTable.ax;
outerHeights = posTable.outerHeights;
outerWidths = posTable.outerWidths;

% Detect number of rows and columns
nRows = round(1 ./ mean(outerHeights));
nColumns = round(1 ./ mean(outerWidths));

% Create a temporary figure with the same number of subplot rows and columns
[figTemp, axTemp] = ...
    create_subplots(nRows, nColumns, 'AlwaysNew', true, ...
                    'Visible', 'off', otherArguments);

% Get the axes positions
axPositionsNew = get(axTemp, 'Position');

% Close the temporary figure
close(figTemp);

% Update axes positions
if nAx > 1
    for iAx = 1:nAx
        set(ax(iAx), 'Position', axPositionsNew{iAx});
    end
else
    set(ax, 'Position', axPositionsNew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%