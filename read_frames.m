function [frames, vidObj] = read_frames (videoPathOrObj, varargin)
%% Reads all frames from a video file
% Usage: [frames, vidObj] = read_frames (videoPathOrObj, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [xylo, v] = read_frames('xylophone.mp4')
%       xyloSome = read_frames('xylophone.mp4', 'TimeWindow', [0.5, 0.9]);
%       xyloHead = read_frames('xylophone.mp4', 'TimeStart', 0.5);
%       xyloTail = read_frames('xylophone.mp4', 'TimeEnd', 0.9);
%       xylo3 = read_frames('xylophone.mp4', 'FrameIndex', 3)
%       xylo3to5 = read_frames('xylophone.mp4', 'IndexStart', 3, 'IndexEnd', 5)
%       xylo60 = read_frames('xylophone.mp4', 'FrameIndex', 60);
%       figure(1); image(xylo60.cdata);
%       xyloAll = read_frames('xylophone.mp4');
%       figure(2); image(xyloAll(60).cdata);
%
% Outputs:
%       frames      - MATLAB movie frame structures, with fields:
%                       cdata    - RGB intensity data
%                       colormap - color map used
%                       time     - relative frame time in seconds
%                       duration - duration of frame in seconds
%                   specified as a structure array
%       vidObj      - VideoReader object
%                   specified as a VideoReader object
%
% Arguments:
%       videoPathOrObj  - path to video file or the VideoReader object
%                       must be a string scalar or a character vector or
%                           a VideoReader object
%       varargin    - 'FrameIndex': frame index to read relative to
%                               time start of the current videoReader object
%                   must be empty or a numeric vector
%                   default == read all frames
%                   - 'IndexStart': first frame index to read relative to
%                               time start of the current videoReader object
%                   must be empty or a numeric vector
%                   default == 1
%                   - 'IndexEnd': last frame index to read relative to
%                               time start of the current videoReader object
%                   must be empty or a numeric vector
%                   default == count_frames(vidObj)
%                   - 'TimeWindow': time window (in sec) to read the frame
%                   must be a 2-element numeric vector
%                   default == entire length
%                   - 'TimeStart': start time (in sec) to read the frame
%                   must be a numeric scalar
%                   default == vidObj.CurrentTime
%                   - 'TimeEnd': end time (in sec) to read the frame
%                   must be a numeric scalar
%                   default == vidObj.Duration
%                   - Any other parameter-value pair for VideoReader()
%
% Requires:
%       cd/create_empty_frames.m
%       cd/create_error_for_nargin.m
%       cd/decide_on_video_object.m
%       cd/struct2arglist.m
%
% Used by:

% File History:
% 2019-09-04 Adapted from https://www.mathworks.com/help/matlab/import_export/read-video-files.html
% 2019-09-05 Added times and durations in the frame structure

%% Hard-coded parameters

%% Default values for optional arguments
frameIndexDefault = [];         % set later
indexStartDefault = [];         % set later
indexEndDefault = [];           % set later
timeWindowDefault = [];         % set later
timeStartDefault = [];          % set later
timeEndDefault = [];            % set later

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
addParameter(iP, 'FrameIndex', frameIndexDefault, ...
    @(x) assert(isempty(x) || isnum(x), ...
                'FrameIndex must be either empty or a numeric scalar!'));
addParameter(iP, 'IndexStart', indexStartDefault, ...
    @(x) assert(isnumericvector(x), ...
                'IndexStart must be either empty or a numeric vector!'));
addParameter(iP, 'IndexEnd', indexEndDefault, ...
    @(x) assert(isnumericvector(x), ...
                'IndexEnd must be either empty or a numeric vector!'));
addParameter(iP, 'TimeWindow', timeWindowDefault, ...
    @(x) assert(isnumericvector(x), ...
                'TimeWindow must be either empty or a numeric vector!'));
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) assert(isnumericvector(x), ...
                'TimeStart must be either empty or a numeric vector!'));
addParameter(iP, 'TimeEnd', timeEndDefault, ...
    @(x) assert(isnumericvector(x), ...
                'TimeEnd must be either empty or a numeric vector!'));

% Read from the Input Parser
parse(iP, videoPathOrObj, varargin{:});
frameIndex = iP.Results.FrameIndex;
indexStart = iP.Results.IndexStart;
indexEnd = iP.Results.IndexEnd;
timeWindow = iP.Results.TimeWindow;
timeStart = iP.Results.TimeStart;
timeEnd = iP.Results.TimeEnd;

% Keep unmatched arguments for the VideoReader() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the video object
vidObj = decide_on_video_object(videoPathOrObj, otherArguments{:});

% Compute the minimum starting time
minTimeStart = 0;

% Compute the maximum end time
maxTimeEnd = vidObj.Duration;

% Decide on the starting index
if isempty(indexStart)
    if ~isempty(frameIndex)
        indexStart = frameIndex;
    else    
        indexStart = 1;
    end
end

% Decide on the ending index
if isempty(indexEnd)
    if ~isempty(frameIndex)
        indexEnd = frameIndex;
    else    
        indexEnd = Inf;
    end
end

% Decide on the starting time
if isempty(timeStart)
    if ~isempty(timeWindow)
        timeStart = timeWindow(1);
    else
        timeStart = minTimeStart;
    end
else
    if ~isempty(timeWindow)
        disp('Warning: TimeStart will override TimeWindow(1)!');
    end
end

% Decide on the end time
if isempty(timeEnd)
    if ~isempty(timeWindow)
        timeEnd = timeWindow(2);
    else
        timeEnd = maxTimeEnd;
    end
else
    if ~isempty(timeWindow)
        disp('Warning: TimeEnd will override TimeWindow(2)!');
    end
end

% Make sure timeStart is valid
if isnan(timeStart) || timeStart < minTimeStart || timeStart >= maxTimeEnd
    timeStart = minTimeStart;
end

% Make sure timeEnd is valid
if isnan(timeEnd) || timeEnd < timeStart || timeEnd >= maxTimeEnd
    timeEnd = maxTimeEnd;
end

%% Do the job
% Read the height and width of the video object
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

% Initialize the MATLAB movie frame structure array
frames = create_empty_frames(vidWidth, vidHeight);

% Initialize the VideoReader object at a specific time
vidObj.CurrentTime = timeStart;

% Initialize a frame index
iFrame = 0;

% Initialize a count for frames to read
count = 0;

% Read frames as long as time is not yet time end and count is not yet index end
while vidObj.CurrentTime < timeEnd && iFrame < indexEnd
    % Increment frame index
    iFrame = iFrame + 1;

    % Read in the frame data if the time is not yet time end
    %   and if the index is within range
    if vidObj.CurrentTime <= timeEnd && ...
            iFrame >= indexStart && iFrame <= indexEnd
        % Increment count for frames to read
        count = count + 1;

        % Read the time of this frame
        %   Note: this must occur before readFrame() is called
        frames(count, 1).time = vidObj.CurrentTime;

        % Read the frame and store in output
        frames(count, 1).cdata = readFrame(vidObj);

        % Read the duration of this frame
        %   Note: this must occur after readFrame() is called
        frames(count, 1).duration = vidObj.CurrentTime - frames(count, 1).time;
    else
        % Read the frame and discard
        readFrame(vidObj);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
