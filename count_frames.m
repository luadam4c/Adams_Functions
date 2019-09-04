function nFrames = count_frames (videoPath, varargin)
%% Count the number of frames in a video file
% Usage: nFrames = count_frames (videoPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       nFrames     - TODO: Description of nFrames
%                   specified as a TODO
%
% Arguments:
%       videoPath     - TODO: Description of videoPath
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/struct2arglist.m
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 201X-XX-XX Created by TODO or Adapted from TODO
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'videoPath');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, videoPath, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
% TODO

%% Preparation
% TODO

%% Do the job
% Create a video object from the video file
vidObj = VideoReader(videoPath);

% Initialize the number of frames
nFrames = 0;

% Increment the frame count
while hasFrame(vidObj)
    % Read a frame and throw it away
    readFrame(vidObj);

    % Increment the frame count
    nFrames = nFrames + 1;
end
    
%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%