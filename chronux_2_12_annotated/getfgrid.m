function [freqs, freqIndices] = getfgrid (fs, nFft, fPass)
%% Helper function to generate frequency grid for spectral estimation
% Usage: [freqs, freqIndices] = getfgrid (fs, nFft, fPass)
% Explanation:
%       Generates the frequency vector and finds indices corresponding 
%       to the requested frequency band (fPass).
%
% Outputs:
%       freqs       - Vector of frequencies
%       freqIndices - Indices in the full grid corresponding to fPass
%
% Arguments:
%       fs          - Sampling frequency
%       nFft        - Number of points in FFT
%       fPass       - Frequency band [fmin fmax]
%
% Requires:
%
% Used by:
%       coherencyc.m
%       coherencycpt.m
%       mtspectrumc.m
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from getfgrid.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 3
    error('Need all arguments');
end

%% Do the job
% Calculate frequency resolution
df = fs / nFft;

% Generate all possible frequencies
allFreqs = 0:df:fs; 

% Truncate to match FFT length (handling potential floating point excess)
allFreqs = allFreqs(1:nFft);

% Determine indices for the passband
if length(fPass) ~= 1
    % Band pass case: [fMin fMax]
    freqIndices = find(allFreqs >= fPass(1) & allFreqs <= fPass(end));
else
    % Single frequency case: find nearest
    [~, freqIndices] = min(abs(allFreqs - fPass));
    % clear fmin (unused)
end

%% Output results
freqs = allFreqs(freqIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%