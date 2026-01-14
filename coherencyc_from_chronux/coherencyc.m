function [coherenceMag, phi, crossSpec, spec1, spec2, freqs, confC, phiStd, coherenceErr] = coherencyc (data1, data2, params)
%% Computes Multi-taper coherency, cross-spectrum and individual spectra
% Usage: [coherenceMag, phi, crossSpec, spec1, spec2, freqs, confC, phiStd, coherenceErr] = coherencyc (data1, data2, params)
% Explanation:
%       Calculates the coherency between two continuous datasets using 
%       multi-taper spectral estimation.
%       Note: Units have to be consistent. See chronux.m for more information.
%
% Outputs:
%       coherenceMag - Magnitude of coherency
%                   (frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       phi          - Phase of coherency
%                   (frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       crossSpec    - Cross spectrum
%                   (frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       spec1        - Spectrum of data1
%                   (frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       spec2        - Spectrum of data2
%                   (frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       freqs        - Frequencies
%       confC        - Confidence level for C at 1-p %
%                   (only for err(1)>=1)
%       phiStd       - Theoretical/Jackknife standard deviation for phi.
%                   (depending on err(1)=1/err(1)=2).
%                   Note that phi + 2 phistd and phi - 2 phistd will give 
%                   95% confidence bands for phi - only for err(1)>=1
%       coherenceErr - Jackknife error bars for C 
%                   (use only for Jackknife - err(1)=2)
%
% Arguments:
%       data1   - First dataset
%                   (in form samples x trials) -- required
%       data2   - Second dataset
%                   (in form samples x trials) -- required
%       params  - Structure with fields:
%           .tapers : precalculated N by K tapers from discrete prolate spheroidal sequences
%                            or in one of the following forms:
%                       (1) A numeric vector [TW K] where TW is the
%                           time-half-bandwidth product and K is the number of
%                           tapers to be used (less than or equal to 2TW-1).
%                       (2) A numeric vector [W T p] where W is the
%                           bandwidth, T is the duration of the data and p
%                           is an integer such that 2TW-p tapers are used.
%                           In this form there is no default i.e. to specify
%                           the bandwidth, you have to specify T and p as well.
%                           Note that the units of W and T have to be consistent:
%                           if W is in Hz, T must be in seconds and vice versa.
%                           Note that these units must also be consistent with
%                           the units of params.Fs: W can be in Hz if and only
%                           if params.Fs is in Hz.
%                       Default: [3 5] (TW = 3, K = 5)
%           .pad    : (padding factor for the FFT) - optional
%                       Can take values -1, 0, 1, 2...
%                       -1 corresponds to no padding.
%                       0 corresponds to padding to the next highest power of 2 etc.
%                       e.g. For N = 500:
%                           if PAD = -1, we do not pad;
%                           if PAD = 0, we pad the FFT to 512 points;
%                           if PAD = 1, we pad to 1024 points etc.
%                       Default: 0
%           .Fs     : (sampling frequency) - optional.
%                       Default: 1
%           .fpass  : (frequency band to be used in the calculation) - optional
%                       In the form [fmin fmax].
%                       Default: all frequencies between 0 and Fs/2
%           .err    : (error calculation) - optional
%                       [1 p] - Theoretical error bars
%                       [2 p] - Jackknife error bars
%                       [0 p] or 0 - no error bars
%                       Default: 0
%           .trialave : (average over trials) - optional
%                       1 : average over trials
%                       0 : don't average
%                       Default: 0
%
% Requires:
%       change_row_to_column.m
%       getparams.m
%       check_consistency.m
%       getfgrid.m
%       dpsschk.m
%       mtfftc.m
%       coherr.m
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from coherencyc.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 2
    error('Need data1 and data2'); 
end

% Ensure column orientation
data1 = change_row_to_column(data1);
data2 = change_row_to_column(data2);

if nargin < 3; params = []; end

% Parse parameters
[tapers, pad, fs, fPass, errConfig, trialAve] = getparams(params);

% Check output argument validity regarding error calculation
if nargout > 8 && errConfig(1) ~= 2 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end
if nargout > 6 && errConfig(1) == 0
    % Errors computed only if err(1) is nonzero. Need to change params and run again.
    error('When errors are desired, err(1) has to be non-zero.');
end

%% Preparation
% Check consistency of data dimensions
nSamples = check_consistency(data1, data2);

% Calculate FFT length
nFft = max(2^(nextpow2(nSamples) + pad), nSamples);

% Get frequency grid
[freqs, freqIndices] = getfgrid(fs, nFft, fPass); 

% Check/Calculate tapers
tapers = dpsschk(tapers, nSamples, fs); 

%% Do the job
% 1. Compute Multi-taper FFTs
j1 = mtfftc(data1, tapers, nFft, fs);
j2 = mtfftc(data2, tapers, nFft, fs);

% 2. Restrict to desired frequencies
j1 = j1(freqIndices, :, :); 
j2 = j2(freqIndices, :, :);

% 3. Compute Spectra and Cross-Spectra
% Average over tapers (dim 2)
crossSpec = squeeze(mean(conj(j1) .* j2, 2));
spec1     = squeeze(mean(conj(j1) .* j1, 2));
spec2     = squeeze(mean(conj(j2) .* j2, 2));

% 4. Handle Trial Averaging
if trialAve
    crossSpec = squeeze(mean(crossSpec, 2)); 
    spec1     = squeeze(mean(spec1, 2)); 
    spec2     = squeeze(mean(spec2, 2)); 
end

% 5. Compute Coherency
c12 = crossSpec ./ sqrt(spec1 .* spec2);
coherenceMag = abs(c12); 
phi = angle(c12);

%% Output results
if nargout >= 9
     [confC, phiStd, coherenceErr] = coherr(coherenceMag, j1, j2, errConfig, trialAve);
elseif nargout == 8
     [confC, phiStd] = coherr(coherenceMag, j1, j2, errConfig, trialAve);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%