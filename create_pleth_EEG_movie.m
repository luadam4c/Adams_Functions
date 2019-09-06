%%function [plotFrames, vidWriter, handles] = create_pleth_EEG_movie (varargin)
%% Creates a synced movie from a .wmv file and a Spike2-exported .mat file in the current directory
% Usage: [plotFrames, vidWriter, handles] = create_pleth_EEG_movie (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_files.m
%       cd/extract_fileparts.m
%       cd/find_in_strings.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/create_time_vectors.m
%       cd/read_frames.m
%       cd/struct2arglist.m
%       cd/write_frames.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-05 Created by Adam Lu
% 

%% Hard-coded parameters

%% TODO: Make optional arguments
spike2MatPath = '';
wmvPath = '';
eegChannelName = 'WIC_2';
% movieType = 'MPEG-4';             % Only works in Windows
movieType = 'Motion JPEG AVI';
outFolder = '';
movieBase = '';

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
% iP = inputParser;
% iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
% parse(iP, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the Spike2-exported mat file
if isempty(spike2MatPath)
    [~, spike2MatPath] = all_files('Extension', 'mat', 'MaxNum', 1, ...
                                    'ForceCellOutput', false);
end

% Decide on the wmv file
[~, wmvPath] = all_files('Extension', 'wmv', 'MaxNum', 1, ...
                                    'ForceCellOutput', false);

% Decide on the output folder
if isempty(outFolder)
    outFolder = extract_fileparts(spike2MatPath, 'directory');
end

% Decide on the movie file base
if isempty(movieBase)
    movieBase = [extract_fileparts(spike2MatPath, 'base'), '_movie'];
end

%% Deal with the Spike2 MATLAB file
% Load .mat file
spike2File = matfile(spike2MatPath);

% Get all the structure names
allStructNames = fieldnames(spike2File);

% Find the structure with EEG trace data
[~, eegStructName] = find_in_strings(eegChannelName, allStructNames);

% Extract the EEG struct
eegStruct = spike2File.(eegStructName);

% Extract the EEG trace data
traceData = eegStruct.values;

% Extract the EEG trace time info
timeStart = eegStruct.start;
siSeconds = eegStruct.interval;
nSamples = eegStruct.length;

% Construct a time vector
tVec = create_time_vectors(nSamples, 'TimeStart', timeStart, ...
                    'SamplingIntervalSeconds', siSeconds, 'TimeUnits', 's', ...
                    'BoundaryMode', 'leftadjust');

%% Deal with the movie file
% Read all frames
[frames, vidObj] = read_frames(wmvPath);

%% Combine into a plot movie
% Create plot movie
[plotFrames, handles] = ...
    create_synced_movie_trace_plot_movie(frames, traceData, 'TimeVec', tVec);

%% Write movie to file
vidWriter = write_frames(plotFrames, 'MovieType', movieType, ...
                'OutFolder', outFolder, 'FileBase', movieBase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%