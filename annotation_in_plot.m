function an = annotation_in_plot (annotationType, xInAxes, yInAxes, varargin)
%% A wrapper function for the annotation() function that accepts x and y values normalized to the axes
% Usage: an = annotation_in_plot (annotationType, xInAxes, yInAxes, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       an          - handle to annotation object
%                   specified as an object handle
% Arguments:
%       annotationType  - line or shape type, see doc for annotation()
%                       must be: see doc for annotation()
%       xInAxes         - [xBegin, xEnd], relative to the axes
%                       must be a 2-element numeric vector
%       yInAxes         - [yBegin yEnd], relative to the axes
%                       must be a 2-element numeric vector
%       varargin    - Any parameter-value pair for the annotation() function
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_numeric.m
%       cd/is_out_of_range.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_pulse_response_with_stimulus.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-28 Added usage of struct2arglist.m
% 2018-12-28 Now deals with positions that are out of range
% 2018-12-28 Now returns empty graphics object if any position is NaN

%% Hard-coded parameters
rangeToCheck = [0, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'annotationType');
addRequired(iP, 'xInAxes', ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addRequired(iP, 'yInAxes', ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, annotationType, xInAxes, yInAxes, varargin{:});

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

% Convert to a cell array
params = struct2arglist(otherArguments);

%% Preparation
% If there is any NaN values, return
if any(isnan(xInAxes) | isnan(yInAxes))
    an = gobjects(1);
    return
end

% Force as column vectors
%   Note: important if flipud() is used later
[xInAxes, yInAxes] = argfun(@force_column_numeric, xInAxes, yInAxes);

% Get the positions of the axes in normalized units relative to the figure
posInFigure = get(gca, 'Position');

% If any position is out of range, annotate partially
if is_out_of_range([xInAxes; yInAxes], rangeToCheck)
    switch annotationType
        case 'line'
            % Replace based on the part that is out of range
            xInAxes(xInAxes < 0) = 0; 
            xInAxes(xInAxes > 1) = 1; 
            yInAxes(yInAxes < 0) = 0; 
            yInAxes(yInAxes > 1) = 1; 
        case {'arrow', 'textarrow'}
            if is_out_of_range([xInAxes(2), yInAxes(2)], rangeToCheck)
                % If the second end is out of range, plot a line instead
                an = annotation_in_plot('line', xInAxes, yInAxes, ...
                                        varargin{:});
                return
            else
                % If the first end is out of range, replace
                xInAxes(xInAxes < 0) = 0; 
                xInAxes(xInAxes > 1) = 1; 
                yInAxes(yInAxes < 0) = 0; 
                yInAxes(yInAxes > 1) = 1; 
            end
        case 'doublearrow'
            if any(xInAxes < 0) && any(xInAxes > 1) || ...
                any(yInAxes < 0) && any(yInAxes > 1)
                % If both ends out of range, draw a line instead
                an = annotation_in_plot('line', xInAxes, yInAxes, ...
                                        varargin{:});
                return
            else
                % If only one end out of range, draw an arrow instead
                % First, reorder x and y values if it is the second end 
                %   that is out of range
                if is_out_of_range([xInAxes(2), yInAxes(2)], rangeToCheck)
                    [xInAxes, yInAxes] = argfun(@flipud, xInAxes, yInAxes);
                end

                % Draw an arrow instead
                an = annotation_in_plot('arrow', xInAxes, yInAxes, ...
                                        varargin{:});
                return
            end
        otherwise
            error('%s is out of range but not implemented yet!!', annotationType);
    end
end

% TODO: If any position is out of range ([0, 1]), expand axes temporarily,
%           draw annotation, fix it, then shift axes back
% if any(xInAxes < 0 | xInAxes > 1 | yInAxes < 0 | yInAxes > 1)
%     % Save old axes
%     xLimitsOld = get(gca, 'XLim');
%     yLimitsOld = get(gca, 'YLim');

%     % Expand axes

%     % Compute the x and y values in normalized units relative to the new axes
% end

% Compute the x and y values in normalized units relative to the figure
xInFigure = posInFigure(1) + posInFigure(3) * xInAxes;
yInFigure = posInFigure(2) + posInFigure(4) * yInAxes;

%% Do the job

% Draw the annotation
an = annotation(annotationType, xInFigure, yInFigure, params{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if xInAxes(1) < 0 || xInAxes(2) > 1 || yInAxes(1) < 0 || yInAxes(2) > 1

if xInAxes(1) < 0 
    xInAxes(1) = 0;
end
if xInAxes(2) > 1 
    xInAxes(2) = 1;
end
if yInAxes(1) < 0 
    yInAxes(1) = 0;
end
if yInAxes(2) > 1
    yInAxes(2) = 1;
end                

if xInAxes(2) < 0 || xInAxes(2) > 1 || ...
        yInAxes(2) < 0 || yInAxes(2) > 1

is_out_of_range (values, rangeToCheck)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%