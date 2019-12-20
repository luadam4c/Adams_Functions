function objects = plot_ball_stick (varargin)
%% Plots a ball-and-stick model
% Usage: objects = plot_ball_stick (radiusSoma, radiusDend, lengthDend, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       objects = plot_ball_stick(3, 1, 10)
%       objects = plot_ball_stick(3, 1, 10, 'FaceColor', 'none')
%       objects = plot_ball_stick(3, 1, 10, 'EdgeColor', 'g')
%
% Outputs:
%       objects     - rectangle objects plotted
%                   specified as rectangle objects
%
% Arguments:
%       radiusSoma  - radius of soma
%                   must be a nonegative scalar
%       radiusDend  - radius of dendrite
%                   must be a nonegative scalar
%       lengthDend  - length of dendrite
%                   must be a nonegative scalar
%       varargin    - 'GeomParams': geometric parameters
%                   must be a scalar struct with fields:
%                       radiusSoma
%                       radiusDendrite
%                       lengthDendrite
%                   default == struct.empty
%                   - 'BallFaceColor': face color for ball
%                   must be a character vector or 3-element numeric array
%                   default == rgb('MediumBlue')
%                   - 'BallEdgeColor': edge color for ball
%                   must be a character vector or 3-element numeric array
%                   default == rgb('MediumBlue')
%                   - 'StickFaceColor': face color for stick
%                   must be a character vector or 3-element numeric array
%                   default == rgb('DarkOrange')
%                   - 'StickEdgeColor': edge color for stick
%                   must be a character vector or 3-element numeric array
%                   default == rgb('DarkOrange')
%                   - 'FaceColor': face color for both ball and stick
%                   must be a character vector or 3-element numeric array
%                   default == [] (not set)
%                   - 'EdgeColor': edge color for both ball and stick
%                   must be a character vector or 3-element numeric array
%                   default == [] (not set)
%                   - Any other parameter-value pair for rectangle()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/hold_on.m
%       cd/hold_off.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/find_passive_params.m
%       cd/m3ha_plot_figure03.m

% File History:
% 2019-12-20 Moved code from find_passive_params.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
radiusSomaDefault = [];
radiusDendDefault = [];
lengthDendDefault = [];
geomParamsDefault = struct.empty;
ballFaceColorDefault = [0, 0, 0.8008];      % rgb('MediumBlue')
ballEdgeColorDefault = [0, 0, 0.8008];      % rgb('MediumBlue')
stickFaceColorDefault = [1, 0.5469, 0];     % rgb('DarkOrange')
stickEdgeColorDefault = [1, 0.5469, 0];     % rgb('DarkOrange')
faceColorDefault = [];
edgeColorDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'radiusSoma', radiusSomaDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addOptional(iP, 'radiusDend', radiusDendDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addOptional(iP, 'lengthDend', lengthDendDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'GeomParams', geomParamsDefault, ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addParameter(iP, 'BallFaceColor', ballFaceColorDefault);
addParameter(iP, 'BallEdgeColor', ballEdgeColorDefault);
addParameter(iP, 'StickFaceColor', stickFaceColorDefault);
addParameter(iP, 'StickEdgeColor', stickEdgeColorDefault);
addParameter(iP, 'FaceColor', faceColorDefault);
addParameter(iP, 'EdgeColor', edgeColorDefault);

% Read from the Input Parser
parse(iP, varargin{:});
radiusSoma = iP.Results.radiusSoma;
radiusDend = iP.Results.radiusDend;
lengthDend = iP.Results.lengthDend;
geomParams = iP.Results.GeomParams;
ballFaceColor = iP.Results.BallFaceColor;
ballEdgeColor = iP.Results.BallEdgeColor;
stickFaceColor = iP.Results.StickFaceColor;
stickEdgeColor = iP.Results.StickEdgeColor;
faceColor = iP.Results.FaceColor;
edgeColor = iP.Results.EdgeColor;

% Keep unmatched arguments for the rectangle() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Replace default colors with user's choice
if ~isempty(faceColor)
    ballFaceColor = faceColor;
    stickFaceColor = faceColor;
end
if ~isempty(edgeColor)
    ballEdgeColor = edgeColor;
    stickEdgeColor = edgeColor;
end

% Read from geomParams if needed
if isempty(radiusSoma) || isempty(radiusDend) || isempty(lengthDend)
    if ~isempty(geomParams)
        radiusSoma = geomParams.radiusSoma;
        radiusDend = geomParams.radiusDendrite;
        lengthDend = geomParams.lengthDendrite;
    else
        error('No geometric parameters passed in!');
    end
end

%% Do the job
% Hold on
wasHold = hold_on;

% Plot soma as a circle
objects(1) = rectangle('Position', radiusSoma * [-1, -1, 2, 2], ...
                        'Curvature', [1, 1], 'FaceColor', ballFaceColor, ...
                        'EdgeColor', ballEdgeColor, ...
                        otherArguments{:});

% Plot dendrite as a rectangle
objects(2) = rectangle('Position', [radiusSoma, -radiusDend, ...
                                    lengthDend, 2 * radiusDend], ...
                        'Curvature', [0, 0], 'FaceColor', stickFaceColor, ...
                        'EdgeColor', stickEdgeColor, ...
                        otherArguments{:});

% Adjust data aspect ratio so that the x and y unit lengths are the same
daspect([1, 1, 1]);

% Hold off
hold_off(wasHold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%