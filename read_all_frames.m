function frames = read_all_frames (videoPathOrObj, varargin)
%% Reads all frames from a video file
% Usage: frames = read_all_frames (videoPathOrObj, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       frames = read_all_frames('xylophone.mp4')
%
% Outputs:
%       frames      - all the frames
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
% 2019-09-04 Adapted from https://www.mathworks.com/help/matlab/import_export/read-video-files.html
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

%% Preparation
% Decide on the video object
vidObj = decide_on_video_object(videoPathOrObj);

%% Do the job
% Read the height and width of the video object
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%