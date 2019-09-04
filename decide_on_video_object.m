function vidObj = decide_on_video_object (videoPathOrObj, varargin)
%% Decide on the video object given path or object
% Usage: vidObj = decide_on_video_object (videoPathOrObj, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       vidObj      - the VideoReader object
%                   specified as a VideoReader object
%
% Arguments:
%       videoPathOrObj  - path to video file or the VideoReader object itself
%                       must be a string scalar or a character vector or
%                           a VideoReader object
%       varargin    - Any other parameter-value pair for VideoReader()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/count_frames.m
%       cd/extract_frame_times.m
%       cd/read_all_frames.m

% File History:
% 2019-09-04 Created by Adam Lu
% 

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

% Read from the Input Parser
parse(iP, videoPathOrObj, varargin{:});

% Keep unmatched arguments for the VideoReader() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
if isa(videoPathOrObj, 'VideoReader')
    % The object is already passed in
    vidObj = videoPathOrObj;
else
    % Create a video object from the video file
    vidObj = VideoReader(videoPathOrObj, otherArguments{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%