function [tapersOutput, eigenValues] = dpsschk (tapersInput, nSamples, fs)
%% Checks or calculates discrete prolate spheroidal sequences tapers ensuring length consistency
% Usage: [tapersOutput, eigenValues] = dpsschk (tapersInput, nSamples, fs)
% Explanation:
%       If precalculated tapers are supplied, checks that they match the 
%       time series length (nSamples). 
%       If parameters [TW K] are supplied, calculates the DPSS tapers.
%
% Outputs:
%       tapersOutput - Calculated or verified tapers (normalized)
%       eigenValues  - Eigenvalues of the tapers
%
% Arguments:
%       tapersInput - Tapers in form:
%                     (i) Precalculated N by K tapers matrix
%                     (ii) [TW K] - time-half-bandwidth product, number of tapers
%       nSamples    - Number of time samples in the data
%       fs          - Sampling frequency (required for normalization)
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from dpsschk.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 3
    error('Need all arguments: tapers, nSamples, fs');
end

%% Do the job
sz = size(tapersInput);

% Check if input is [TW K] (1x2 vector)
if sz(1) == 1 && sz(2) == 2
    % Calculate DPSS tapers
    % tapersInput(1) = TW (Time-Half-Bandwidth product)
    % tapersInput(2) = K (Number of tapers)
    [tapersOutput, eigenValues] = dpss(nSamples, tapersInput(1), tapersInput(2));
    
    % Normalize tapers:
    % dpss computes tapers such that SUM of squares = 1.
    % We need integral of squares = 1, so multiply by sqrt(Fs).
    tapersOutput = tapersOutput * sqrt(fs);

% Check if precalculated tapers match data length
elseif nSamples ~= sz(1)
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
else
    % Tapers are already precalculated and correct size
    tapersOutput = tapersInput;
    eigenValues = []; % No eigenvalues returned for pre-calc check
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%