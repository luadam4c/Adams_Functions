function [spectData, freqHz, timeInstantsSeconds] = ...
                compute_spectrogram (dataValues, siSeconds, varargin)
%% Computes a spectrogram
% Usage: [spectData, freqHz, timeInstantsSeconds] = ...
%               compute_spectrogram (dataValues, siSeconds, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [data, freq, time] = compute_spectrogram(rand(100, 1), 0.1);
%       [data, freq, time] = compute_spectrogram(rand(1000, 1), 0.005);
%       plot_spectrogram(data, time, freq);
%
% Outputs:
%       spectData               - short-time Fourier transform
%                                   in [data units x seconds]
%                               specified as a numeric matrix
%       freqHz                  - cyclical frequencies in Hz
%                               specified as a numeric vector
%       timeInstantsSeconds     - midpoints of each time bin in seconds
%                               specified as a numeric vector
%
% Arguments:
%       dataValues  - data values
%                   must be a numeric vector readable 
%                       by the spectrogram() function
%       siSeconds   - sampling interval in seconds
%                   must be a positive scalar
%       varargin    - 'BinWidthSeconds': bin width in seconds
%                   must be a numeric scalar
%                   default == 1
%                   - 'OverlapSeconds': bin overlap in seconds
%                   must be empty or a numeric scalar
%                   default == half of bin width
%                   - 'StartTimeSeconds': data start time in seconds
%                   must be a nonnegative scalar
%                   default == 0
%                   - Any other parameter-value pair for spectrogram()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/plot_spectrogram_multiunit.m
%       cd/plot_traces_spike2_mat.m

% File History:
% 2019-10-15 Moved from plot_traces_spike2_mat.m
% 2020-08-05 Made 'OverlapSeconds' an optional argument
% 2020-08-05 Made 'StartTimeSeconds' an optional argument

%% Hard-coded parameters

%% Default values for optional arguments
binWidthSecondsDefault = 1;     % 1 second bins by default
overlapSecondsDefault = [];     % set later
startTimeSecondsDefault = 0;    % data starts at 0 seconds by default

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
addRequired(iP, 'dataValues');
addRequired(iP, 'siSeconds', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BinWidthSeconds', binWidthSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'OverlapSeconds', overlapSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'StartTimeSeconds', startTimeSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, dataValues, siSeconds, varargin{:});
binWidthSeconds = iP.Results.BinWidthSeconds;
overlapSeconds = iP.Results.OverlapSeconds;
startTimeSeconds = iP.Results.StartTimeSeconds;

% Keep unmatched arguments for the spectrogram() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Compute the bin width in samples
binWidthSamples = round(binWidthSeconds / siSeconds);

% Set a default overlap in samples
if isempty(overlapSeconds)
    overlapSamples = round(binWidthSamples / 2);
else
    overlapSamples = round(overlapSeconds / siSeconds);
end

% Compute the sampling frequency in Hz
samplingFreqHz = 1 / siSeconds;

%% Do the job
% Compute the spectrogram
%   Note: time instants are the midpoints of each time window
%           the spectrogram values are complex values 
%           in [data units x seconds]
[spectData, freqHz, timeInstantsSecondsRel] = ...
    spectrogram(dataValues, binWidthSamples, ...
                overlapSamples, [], samplingFreqHz, otherArguments{:});

% Add start time to time instants
timeInstantsSeconds = timeInstantsSecondsRel + startTimeSeconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%