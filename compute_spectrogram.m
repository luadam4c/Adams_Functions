function [spectData, freqHz, timeInstantsSeconds] = ...
                compute_spectrogram (eegValues, siSeconds, varargin)
%% Computes a spectrogram
% Usage: [spectData, freqHz, timeInstantsSeconds] = ...
%               compute_spectrogram (eegValues, siSeconds, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [data, freq, time] = compute_spectrogram(rand(100, 1), 0.1);
%
% Outputs:
%       spectData               - TODO: Description of spectData
%                               specified as a TODO
%       freqHz                  - TODO: Description of freqHz
%                               specified as a TODO
%       timeInstantsSeconds     - TODO: Description of timeInstantsSeconds
%                               specified as a TODO
%
% Arguments:
%       eegValues   - TODO: Description of eegValues
%                   must be a TODO
%       siSeconds   - TODO: Description of siSeconds
%                   must be a TODO
%       varargin    - 'BinWidthSeconds': bin width in seconds
%                   must be a numeric scalar
%                   default == 1
%                   - Any other parameter-value pair for spectrogram()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/plot_traces_spike2_mat.m

% File History:
% 2019-10-15 Moved from plot_traces_spike2_mat.m
% 

%% Hard-coded parameters
% TODO: Make optional arguments
overlapSeconds = [];

%% Default values for optional arguments
binWidthSecondsDefault = 1;     % 1 second bins by default

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
addRequired(iP, 'eegValues');
addRequired(iP, 'siSeconds');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BinWidthSeconds', binWidthSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, eegValues, siSeconds, varargin{:});
binWidthSeconds = iP.Results.BinWidthSeconds;

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
[spectData, freqHz, timeInstantsSeconds] = ...
    spectrogram(eegValues, binWidthSamples, ...
                overlapSamples, [], samplingFreqHz, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%