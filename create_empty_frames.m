function frames = create_empty_frames (width, height, varargin)
%% Creates an empty MATLAB movie frames
% Usage: frames = create_empty_frames (width, height, dimensions (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       create_empty_frames(3, 2)
%       create_empty_frames(3, 2, [5, 1])
%       create_empty_frames(3, 2, [10, 2])
%       create_empty_frames(3, 2, [10, 2], 'Duration', 1)
%
% Outputs:
%       frames     	- empty MATLAB movie frames structure, with fields:
%                       cdata    - RGB intensity data
%                       colormap - color map used
%                       time     - time of the frames in seconds
%                       duration - duration of frames in seconds
%                   specified as a structure array
%
% Arguments:
%       width       - width of the frames
%                   must be a positive scalar
%       height      - height of the frames
%                   must be a positive scalar
%       dimensions  - (opt) dimensions for an array
%                   must be a positive integer 2-element vector
%       varargin    - 'Duration': frame duration
%                   must be a positive scalar
%                   default == NaN
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/read_frames.m

% File History:
% 2019-09-04 Adapted from https://www.mathworks.com/help/matlab/import_export/read-video-files.html
% 2019-09-05 Added times and durations in the frames structure
% 2019-09-05 Added dimensions as an optional argument
% 2019-09-05 Added 'Duration' as an optional argument
% TODO: Add 'ColorMap' as an optional argument

%% Default values for optional arguments
dimensionsDefault = [];
durationDefault = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'width', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addRequired(iP, 'height', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Add optional arguments to the Input Parser
addOptional(iP, 'dimensions', dimensionsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Duration', durationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, height, width, varargin{:});
dimensions = iP.Results.dimensions;
duration = iP.Results.Duration;

%% Do the job
% Create one empty frame
frames = struct('cdata', zeros(height, width, 3, 'uint8'), ...
                'colormap', [], 'time', NaN, 'duration', duration);

% Expand to an array
if ~isempty(dimensions)
    frames = repmat(frames, dimensions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
