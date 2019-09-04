function nFrames = count_frames (videoPathOrObj, varargin)
%% Count the number of frames in a video file
% Usage: nFrames = count_frames (videoPathOrObj, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       count_frames('xylophone.mp4')
%
% Outputs:
%       nFrames     - number of frames in the video
%                   specified as a positive integer
%
% Arguments:
%       videoPathOrObj  - path to video file or the VideoReader object
%                       must be a string scalar or a character vector or
%                           a VideoReader object
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for VideoReader()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/decide_on_video_object.m
%       cd/struct2arglist.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-03 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'videoPathOrObj');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, videoPathOrObj, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the VideoReader() function
otherArguments = struct2arglist(iP.Unmatched);

%% Prepration
% Decide on the video object
vidObj = decide_on_video_object(videoPathOrObj, otherArguments{:});

%% Do the job
% Initialize the number of frames
nFrames = 0;

% Increment the frame count
while hasFrame(vidObj)
    % Read a frame and throw it away
    readFrame(vidObj);

    % Increment the frame count
    nFrames = nFrames + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%