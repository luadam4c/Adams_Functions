function jFft = mtfftc (data, tapers, nFft, fs)
%% Multi-taper fourier transform for continuous data
% Usage: jFft = mtfftc (data, tapers, nFft, fs)
% Explanation:
%       Computes the FFT of data multiplied by Slepian tapers.
%
% Outputs:
%       jFft    - FFT in form: frequency index x taper index x channels/trials
%
% Arguments:
%       data    - samples x channels/trials (or single vector)
%       tapers  - Precalculated tapers
%       nFft    - Length of padded data
%       fs      - Sampling frequency
%
% Requires:
%       change_row_to_column.m
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from mtfftc.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 4
    error('Need all input arguments');
end

%% Preparation
% Ensure data is column-oriented
data = change_row_to_column(data);

% Get dimensions
[nSamples, nChannels] = size(data);
[nTapers, nK] = size(tapers); % nTapers should match nSamples

% Check compatibility
if nTapers ~= nSamples
    error('length of tapers is incompatible with length of data');
end

%% Do the job
% Replicate tapers for each channel: [Samples, Tapers, Channels]
% tapers currently [Samples, Tapers]
tapersObs = tapers(:, :, ones(1, nChannels)); 

% Replicate data for each taper: [Samples, Tapers, Channels]
% data currently [Samples, Channels] -> add Taper dimension
dataObs = data(:, :, ones(1, nK)); 

% Permute data to match tapers dimensions for element-wise multiplication
% Current dataObs: [Samples, Channels, Tapers]
% Desired: [Samples, Tapers, Channels]
dataObs = permute(dataObs, [1 3 2]); 

% Multiply data by tapers
dataProj = dataObs .* tapersObs; 

% Compute FFT
% Result is normalized by sampling frequency
jFft = fft(dataProj, nFft) / fs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%