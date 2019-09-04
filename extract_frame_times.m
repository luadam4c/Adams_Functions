function frameTimes = extract_frame_times (videoPathOrObj, varargin)
%% Extracts all the frame start times in a video file
% Usage: frameTimes = extract_frame_times (videoPathOrObj, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       extract_frame_times('xylophone.mp4')
%
% Outputs:
%       frameTimes  - all the frame times
%                   specified as a numeric vector
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
% 2019-09-04 Created by Adam Lu
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
% Initialize the VideoReader object at time zero
vidObj.CurrentTime = 0;

% Initialize a frame count
count = 0;

% Initialize frame times
frameTimes = [];

% Record all times
while hasFrame(vidObj)
    % Increment the frame count
    count = count + 1;

    % Read the time of this frame
    frameTimes(count, 1) = vidObj.CurrentTime;

    % Read this frame and throw it away
    readFrame(vidObj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%