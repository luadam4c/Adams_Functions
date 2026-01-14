function [coherenceMag, phi, crossSpec, spec1, spec2, freqs, zeroSp, confC, phiStd, coherenceErr] = coherencycpt (data1, data2, params, fsCorr, timeGrid)
%% Multi-taper coherency, cross-spectrum and individual spectra (Continuous vs Point Process)
% Usage: [coherenceMag, phi, crossSpec, spec1, spec2, freqs, zeroSp, confC, phiStd, coherenceErr] = coherencycpt (data1, data2, params, fsCorr, timeGrid)
% Explanation:
%       Calculates the coherency between a continuous dataset (data1) and a 
%       point process dataset (data2, e.g., spike times) using multi-taper 
%       spectral estimation.
%
% Outputs:
%       coherenceMag - Magnitude of coherency
%                   (frequencies x trials if trialAve=0; dimension frequencies if trialAve=1)
%       phi          - Phase of coherency
%                   (frequencies x trials if trialAve=0; dimension frequencies if trialAve=1)
%       crossSpec    - Cross spectrum
%                   (frequencies x trials if trialAve=0; dimension frequencies if trialAve=1)
%       spec1        - Spectrum of data1
%                   (frequencies x trials if trialAve=0; dimension frequencies if trialAve=1)
%       spec2        - Spectrum of data2
%                   (frequencies x trials if trialAve=0; dimension frequencies if trialAve=1)
%       freqs        - Frequencies
%       zeroSp       - Indicator for trials where no spikes were found
%                   (1 for zero spikes, 0 otherwise)
%       confC        - Confidence level for C at 1-p %
%                   (only for err(1)>=1)
%       phiStd       - Theoretical/Jackknife standard deviation for phi.
%                   (depending on err(1)=1/err(1)=2).
%       coherenceErr - Jackknife error bars for C 
%                   (use only for Jackknife - err(1)=2)
%
% Arguments:
%       data1       - Continuous data
%                   (in form samples x trials) -- required
%       data2       - Point process data (Spikes)
%                   (structure array of spike times with dimension trials; 
%                   also accepts 1d array of spike times) -- required
%       params      - Structure with fields tapers, pad, Fs, fpass, err, trialave
%                   (see function_template or coherencyc for details) -- optional
%       fsCorr      - Finite size corrections
%                   0 (don't use) or 1 (use).
%                   (available only for spikes). Default: 0
%       timeGrid    - Time grid over which the tapers are to be calculated
%                   (Useful when calling from a moving window routine). 
%                   If empty, spike times are used to define the grid.
%
% Requires:
%       check_consistency.m
%       coherr.m
%       dpsschk.m
%       getfgrid.m
%       getparams.m
%       mtfftc.m
%       mtfftpt.m
%
% Used by:
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from coherencycpt.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 2
    error('Need data1 and data2'); 
end
if nargin < 3; params = []; end
if nargin < 4 || isempty(fsCorr); fsCorr = 0; end

[tapers, pad, fs, fPass, errConfig, trialAve, params] = getparams(params);

% Determine time grid if not provided
if nargin < 5 || isempty(timeGrid)
    [nPoints, ~] = size(data1);
    dt = 1/fs;
    timeGrid = 0:dt:(nPoints-1)*dt; % time grid for prolates
end

% Check output argument validity regarding error calculation
if nargout > 7 && errConfig(1) == 0
    error('When errors are desired, err(1) has to be non-zero.');
end
if nargout > 9 && errConfig(1) ~= 2
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end

%% Preparation
% Check consistency
[~, nChannels] = check_consistency(data1, data2, 1);

% Initialize zero-spike indicator
zeroSp = zeros(1, nChannels); 

% Define lengths for FFT
nPointsForDpss = length(timeGrid); 
nFft = max(2^(nextpow2(nPointsForDpss) + pad), nPointsForDpss); 

% Get frequency grid and tapers
[freqs, freqIndices] = getfgrid(fs, nFft, fPass); 
tapers = dpsschk(tapers, nPointsForDpss, fs); 

%% Do the job
% 1. Fourier Transform of Continuous Data (data1)
j1 = mtfftc(data1, tapers, nFft, fs); 
j1 = j1(freqIndices, :, :); % restrict to required frequencies

% 2. Fourier Transform of Discrete Data (data2)
[j2, ~, nSpikes2] = mtfftpt(data2, tapers, nFft, timeGrid, freqs, freqIndices); 

% 3. Identify empty trials
zeroSp(nSpikes2 == 0) = 1; 

% 4. Compute Spectra and Cross-Spectra
% Average over tapers (dim 2)
crossSpec = squeeze(mean(conj(j1) .* j2, 2)); 
spec1     = squeeze(mean(conj(j1) .* j1, 2)); 
spec2     = squeeze(mean(conj(j2) .* j2, 2)); 

% 5. Handle Trial Averaging
if trialAve
    crossSpec = squeeze(mean(crossSpec, 2)); 
    spec1     = squeeze(mean(spec1, 2)); 
    spec2     = squeeze(mean(spec2, 2)); 
end

% 6. Compute Coherency
coherencyComplex = crossSpec ./ sqrt(spec1 .* spec2);
coherenceMag = abs(coherencyComplex);
phi = angle(coherencyComplex);

%% Output results
if nargout >= 9
    if fsCorr == 1
        [confC, phiStd, coherenceErr] = coherr(coherenceMag, j1, j2, errConfig, trialAve, [], nSpikes2); 
    else
        [confC, phiStd, coherenceErr] = coherr(coherenceMag, j1, j2, errConfig, trialAve); 
    end
elseif nargout == 8
    if fsCorr == 1
        [confC, phiStd] = coherr(coherenceMag, j1, j2, errConfig, trialAve, [], nSpikes2); 
    else
        [confC, phiStd] = coherr(coherenceMag, j1, j2, errConfig, trialAve); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%