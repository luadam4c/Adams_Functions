function writer = write_frames (frames, varargin)
%% Write frames to a file
% Usage: writer = write_frames (frames, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [xylo, v] = read_frames('xylophone.mp4')
%       write_frames(xylo, 'VideoObject', v)
%
% Outputs:
%       writer      - the VideoWriter object created for the video
%                   specified as a VideoWriter object
%
% Arguments:
%       frames      - MATLAB movie frame structures, with fields:
%                       cdata    - RGB intensity data
%                       colormap - color map used
%                       time     - relative frame time in seconds
%                       duration - duration of frame in seconds
%                   must be a structure array
%       varargin    - 'FrameRate': frame rate in Hz
%                   must be a positive scalar
%                   default == detected from frames or VideoObject, or 12 Hz
%                   - 'MovieType': movie type
%                   must be an unambiguous, case-insensitive match to one of: 
%                     'Archival'          - Motion JPEG 2000 file 
%                                             with lossless compression
%                     'Motion JPEG AVI'   - AVI file using Motion JPEG encoding
%                     'Motion JPEG 2000'  - Motion JPEG 2000 file
%                     'MPEG-4'            - MPEG-4 file with H.264 encoding 
%                                             (systems with Windows 7 or later, 
%                                             or macOS 10.7 and later)
%                     'Uncompressed AVI'  - Uncompressed with RGB24 video
%                     'Indexed AVI'       - Uncompressed with indexed video
%                     'Grayscale AVI'     - Uncompressed with grayscale video
%                   default == 'Motion JPEG AVI'
%                   - 'VideoObject': VideoReader object for the video
%                   must be a VideoReader object
%                   default == VideoReader.empty
%                   - 'OutFolder': directory to place .atf file
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileBase': file base for the movie
%                   must be a string scalar or a character vector
%                   default == 'newVideo'
%                   - Any other parameter-value pair for VideoWriter()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_pleth_EEG_movie.m

% File History:
% 2019-09-05 Created by Adam Lu
% TODO: Make OutFolder, FileBase, MovieType optional arguments
% 

%% Hard-coded parameters
validMovieTypes = {'Archival', 'Motion JPEG AVI', 'Motion JPEG 2000', ...
                    'MPEG-4', 'Uncompressed AVI', 'Indexed AVI', ...
                    'Grayscale AVI'};

% TODO: Make optional arguments
extraFields = {'time', 'duration'};

%% Default values for optional arguments
frameRateDefault = [];
movieTypeDefault = 'Motion JPEG AVI';
videoObjectDefault = VideoReader.empty;
outFolderDefault = pwd;
fileBaseDefault = 'newVideo';

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
addRequired(iP, 'frames');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FrameRate', frameRateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'MovieType', movieTypeDefault);
    @(x) any(validatestring(x, validMovieTypes));
addParameter(iP, 'VideoObject', videoObjectDefault, ...
    @(x) validateattributes(x, {'VideoReader'}, {'2d'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, frames, varargin{:});
frameRate = iP.Results.FrameRate;
movieType = validatestring(iP.Results.MovieType, validMovieTypes);
vidObj = iP.Results.VideoObject;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;

% Keep unmatched arguments for the VideoWriter() function
otherArguments = struct2arglist(iP.Unmatched);

% Decide on the frame rate in Hz
if isempty(frameRate)
    if isfield(frames(1), 'duration')
        frameRate = 1 / frames(1).duration;
    elseif ~isempty(vidObj)
        frameRate = vidObj.FrameRate;
    else
        frameRate = 12;
    end
end

% Remove 'time' and 'duration' to match MATLAB's frames structure
framesMatlab = rmfield_custom(frames, extraFields);
    
% Create a path for the movie
moviePathBase = fullfile(outFolder, fileBase);

% Create a VideoWriter object
writer = VideoWriter(moviePathBase, movieType, otherArguments{:});

% Set the frame rate in Hz
writer.FrameRate = frameRate;

% Open the VideoWriter object
open(writer);

% Write frames to the file
writeVideo(writer, framesMatlab);

% Close the VideoWriter object
close(writer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outStruct = rmfield_custom (inStruct, fieldNames)
%% Removes a field from a structure only if it exists
% TODO: Pull out as rmfield_custom.m

% Initialize as the input structure
outStruct = inStruct;

% Remove fields one at a time
if iscell(fieldNames)
    for iField = 1:numel(fieldNames)
        % Remove the field if it exists
        outStruct = rmfield_if_exists(outStruct, fieldNames{iField});
    end
else
    outStruct = rmfield_if_exists(outStruct, fieldNames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outStruct = rmfield_if_exists (inStruct, fieldName)
%% Removes a field from a structure only if it exists
% TODO: Pull out as part of rmfield_custom.m

if isfield(inStruct, fieldName)
    outStruct = rmfield(inStruct, fieldName);
else
    outStruct = inStruct;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%