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
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_pulse_response_with_stimulus.m

% File History:
% 2018-12-19 Created by Adam Lu
% 

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

%% Preparation
% Get the positions of the axes in normalized units relative to the figure
posInFigure = get(gca, 'Position');

% Compute the x and y values in normalized units relative to the figure
xInFigure = posInFigure(1) + posInFigure(3) * xInAxes;
yInFigure = posInFigure(2) + posInFigure(4) * yInAxes;

%% Do the job

% TODO: make a function struct2arglist.m
names = fieldnames(otherArguments);

values = struct2cell(otherArguments);

params = force_column_cell(transpose([names, values]), 'ToLinearize', true);

% Draw the annotation
an = annotation(annotationType, xInFigure, yInFigure, params{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%