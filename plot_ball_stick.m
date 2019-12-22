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
%                   - 'BallColor': ball color for both edge and face
%                   must be a character vector or 3-element numeric array
%                   default == [] (not set)
%                   - 'StickColor': stick color for both edge and face
%                   must be a character vector or 3-element numeric array
%                   default == [] (not set)
%                   - 'FaceColor': face color for both ball and stick
%                   must be a character vector or 3-element numeric array
%                   default == [] (not set)
%                   - 'EdgeColor': edge color for both ball and stick
%                   must be a character vector or 3-element numeric array
%                   default == [] (not set)
%                   - 'BallCurvature': curvature for the ball
%                   must be a character vector or 3-element numeric array
%                   default == [1, 1] (circle)
%                   - 'NStickSegments': number of stick segments
%                   must be a positive integer scalar
%                   default == 1
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
% 2019-12-21 Added 'BallCurvature' as an optional argument
% 2019-12-21 Added 'NStickSegments' as an optional argument

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
ballColorDefault = [];
stickColorDefault = [];
faceColorDefault = [];
edgeColorDefault = [];
ballCurvatureDefault = [1, 1];
nStickSegmentsDefault = 1;

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
addParameter(iP, 'BallColor', ballColorDefault);
addParameter(iP, 'StickColor', stickColorDefault);
addParameter(iP, 'FaceColor', faceColorDefault);
addParameter(iP, 'EdgeColor', edgeColorDefault);
addParameter(iP, 'BallCurvature', ballCurvatureDefault);
addParameter(iP, 'NStickSegments', nStickSegmentsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

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
ballColor = iP.Results.BallColor;
stickColor = iP.Results.StickColor;
faceColor = iP.Results.FaceColor;
edgeColor = iP.Results.EdgeColor;
ballCurvature = iP.Results.BallCurvature;
nStickSegments = iP.Results.NStickSegments;

% Keep unmatched arguments for the rectangle() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Replace ball or stick colors with user's choice
if ~isempty(ballColor)
    ballFaceColor = ballColor;
    ballEdgeColor = ballColor;
end
if ~isempty(stickColor)
    stickFaceColor = stickColor;
    stickEdgeColor = stickColor;
end

% Replace face or edge colors with user's choice
%   Note: This has precedence over ballColor or stickColor
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

% Compute the length of each dendritic segment
lengthSegment = lengthDend / nStickSegments;

%% Do the job
% Hold on
wasHold = hold_on;

% Plot soma as a circle
objects(1) = rectangle('Position', radiusSoma * [-1, -1, 2, 2], ...
                    'Curvature', ballCurvature, 'FaceColor', ballFaceColor, ...
                    'EdgeColor', ballEdgeColor, ...
                    otherArguments{:});

% Plot dendrite(s) as a rectangle
for iSegment = 1:nStickSegments
    % Compute the left end point for the dendritic segment
    posLeft = radiusSoma + lengthSegment * (iSegment - 1);

    % Plot the dendritic segment as a rectangle
    objects(1 + iSegment) = ...
        rectangle('Position', [posLeft, -radiusDend, ...
                                lengthSegment, 2 * radiusDend], ...
                    'Curvature', [0, 0], 'FaceColor', stickFaceColor, ...
                    'EdgeColor', stickEdgeColor, ...
                    otherArguments{:});
end

% Adjust data aspect ratio so that the x and y unit lengths are the same
daspect([1, 1, 1]);

% Hold off
hold_off(wasHold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%