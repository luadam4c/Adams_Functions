function [spectrum, freqs, specErr] = mtspectrumc (data, params)
%% Multi-taper spectrum - continuous process
% Usage: [spectrum, freqs, specErr] = mtspectrumc (data, params)
% Explanation:
%       Calculates the spectrum of a continuous dataset using multi-taper
%       estimation.
%
% Outputs:
%       spectrum    - Spectrum 
%                   (frequency x channels/trials if trialAve=0; 
%                    frequency if trialAve=1)
%       freqs       - Frequencies
%       specErr     - Error bars (only for err(1)>=1)
%
% Arguments:
%       data        - Continuous data
%                   (in form samples x channels/trials) -- required
%       params      - Structure with fields tapers, pad, Fs, fpass, err, trialave
%                   -- optional
%
% Requires:
%       change_row_to_column.m
%       dpsschk.m
%       getfgrid.m
%       getparams.m
%       mtfftc.m
%       specerr.m
%
% Used by:
%
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from mtspectrumc.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 1
    error('Need data'); 
end
if nargin < 2; params = []; end

[tapers, pad, fs, fPass, errConfig, trialAve, params] = getparams(params);

if nargout > 2 && errConfig(1) == 0 
    error('When specErr is desired, err(1) has to be non-zero.');
end

% Ensure column orientation
data = change_row_to_column(data);

%% Preparation
nSamples = size(data, 1);
nFft = max(2^(nextpow2(nSamples) + pad), nSamples);

% Get frequency grid
[freqs, freqIndices] = getfgrid(fs, nFft, fPass); 

% Check tapers
tapers = dpsschk(tapers, nSamples, fs); 

%% Do the job
% Compute Multi-taper FFT
jData = mtfftc(data, tapers, nFft, fs);

% Restrict to desired frequencies
jData = jData(freqIndices, :, :);

% Compute Spectrum (Mean of conjugate product across tapers)
% Permute to move frequency to first dimension
spectrum = permute(mean(conj(jData) .* jData, 2), [1 3 2]);

% Handle trial averaging
if trialAve
    spectrum = squeeze(mean(spectrum, 2));
else
    spectrum = squeeze(spectrum);
end

%% Output results
if nargout == 3 
   specErr = specerr(spectrum, jData, errConfig, trialAve);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%