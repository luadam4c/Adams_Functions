function [plotFrames, vidWriter, handles] = create_pleth_EEG_movies (varargin)
%% Creates a synced movie from a .wmv file and a Spike2-exported .mat file in the current directory
% Usage: [plotFrames, vidWriter, handles] = create_pleth_EEG_movies (varargin)
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
%       cd/argfun.m
%       cd/all_files.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/create_time_vectors.m
%       cd/parse_spike2_mat.m
%       cd/read_frames.m
%       cd/struct2arglist.m
%       cd/write_frames.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-05 Created by Adam Lu
% 2019-09-06 Added pleth and EMG channels
% 2019-10-15 Now uses find_matching_files.m
% 2019-10-15 Now uses parse_spike2_mat.m
% TODO: Add Light on info as raster
% TODO: Add spectrogram
% TODO: Make optional arguments

%% Hard-coded parameters

%% TODO: Make optional arguments
directory = pwd;
spike2MatPaths = {};
wmvPaths = {};
channelNames = {'Pleth 2'; 'WIC#2'; 'WIC#1'};
traceLabels = {'Pleth (mL/sec)'; 'EEG amp (uV)'; 'EMG amp (uV)'};
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
if isempty(spike2MatPaths)
    [~, spike2MatPaths] = ...
        all_files('Directory', directory, 'Extension', 'mat', ...
                    'ForceCellOutput', true);
end

% Decide on the wmv file
if isempty(wmvPaths)
    [~, wmvPaths] = ...
        find_matching_files(spike2MatPaths, 'PartType', 'keyword', ...
                            'Extension', 'wmv', 'ForceCellOutput', true);
end

%% Do the job
[plotFrames, vidWriter, handles] = ...
    cellfun(@(x, y) create_one_pleth_EEG_movie(x, y, channelNames, ...
                            traceLabels, outFolder, movieBase, movieType), ...
            spike2MatPaths, wmvPaths, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [plotFrames, vidWriter, handles] = ...
                create_one_pleth_EEG_movie (spike2MatPath, wmvPath, ...
                                        channelNames, traceLabels, ...
                                        outFolder, movieBase, movieType)

% Decide on the output folder
if isempty(outFolder)
    outFolder = extract_fileparts(spike2MatPath, 'directory');
end

% Decide on the movie file base
if isempty(movieBase)
    movieBase = [extract_fileparts(spike2MatPath, 'base'), '_movie'];
end

%% Deal with the Spike2 MATLAB file
% Load channel data
spike2Table = parse_spike2_mat(spike2MatPath, 'ChannelNames', channelNames);

% Extract the trace data
traceData = spike2Table.channelValues;

% Extract the trace time info
timeStart = spike2Table.channelStarts;
siSeconds = spike2Table.siSeconds;
nSamples = spike2Table.nSamples;

% Compute the average
[timeStart, siSeconds, nSamples] = ...
    argfun(@nanmean, timeStart, siSeconds, nSamples);

% Construct a time vector
tVec = create_time_vectors(nSamples, 'TimeStart', timeStart, ...
                    'SamplingIntervalSeconds', siSeconds, 'TimeUnits', 's', ...
                    'BoundaryMode', 'leftadjust');

%% Deal with the movie file
% Read all frames
frames = read_frames(wmvPath);

%% Combine into a plot movie
% Create plot movie
[plotFrames, handles] = ...
    create_synced_movie_trace_plot_movie(frames, traceData, 'TimeVec', tVec, ...
                                            'TraceLabels', traceLabels);

%% Write movie to file
vidWriter = write_frames(plotFrames, 'MovieType', movieType, ...
                'OutFolder', outFolder, 'FileBase', movieBase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Decide on the Spike2-exported mat file
if isempty(spike2MatPath)
    [~, spike2MatPath] = all_files('Extension', 'mat', 'MaxNum', 1, ...
                                    'ForceCellOutput', false);
end

% Decide on the wmv file
% TODO: Use find_matching_files
if isempty(wmvPath)
    [~, wmvPath] = all_files('Extension', 'wmv', 'MaxNum', 1, ...
                                        'ForceCellOutput', false);
end

eegChannelName = 'WIC_2';
emgChannelName = 'WIC_1';
plethChannelName = 'Pleth_2';

% Load .mat file
spike2File = matfile(spike2MatPath);

% Get all the structure names
allStructNames = fieldnames(spike2File);

% Find the structure with trace data
[~, plethStructName] = find_in_strings(plethChannelName, allStructNames);
[~, eegStructName] = find_in_strings(eegChannelName, allStructNames);
[~, emgStructName] = find_in_strings(emgChannelName, allStructNames);

% Extract the structures
plethStruct = spike2File.(plethStructName);
eegStruct = spike2File.(eegStructName);
emgStruct = spike2File.(emgStructName);

% Extract the trace data
traceData = cell(3, 1);
traceData{1} = plethStruct.values;
traceData{2} = eegStruct.values;
traceData{3} = emgStruct.values;

% Extract the trace time info
timeStart = eegStruct.start;
siSeconds = eegStruct.interval;
nSamples = eegStruct.length;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
