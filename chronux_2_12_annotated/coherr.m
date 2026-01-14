function [confC, phiStd, coherenceErr] = coherr (coherenceMag, j1, j2, errConfig, trialAve, numSpikes1, numSpikes2)
%% Computes confidence intervals for coherency
% Usage: [confC, phiStd, coherenceErr] = coherr (coherenceMag, j1, j2, errConfig, trialAve, numSpikes1, numSpikes2)
% Explanation:
%       Computes lower and upper confidence intervals on the coherency 
%       given the tapered fourier transforms.
%       Supports asymptotic estimates and Jackknife estimates.
%
% Outputs:
%       confC        - Confidence level for Coherence (only if errConfig(1) >= 1)
%       phiStd       - Theoretical or Jackknife standard deviation for phase
%       coherenceErr - Jackknife error bars for Coherence (only if errConfig(1) == 2)
%
% Arguments:
%       coherenceMag - Coherence magnitude
%       j1, j2       - Tapered fourier transforms
%       errConfig    - [errType p] 
%                      (errType=1: asymptotic, errType=2: Jackknife)
%                      (p: p-value for error estimates)
%       trialAve     - 0: no averaging, 1: perform trial averaging
%       numSpikes1   - (opt) number of spikes for data1 (finite size corrections)
%       numSpikes2   - (opt) number of spikes for data2
%
% Requires:
%
% Used by:
%       coherencyc.m
%       coherencycpt.m
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from coherr.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 5
    error('Need at least 5 input arguments'); 
end
if errConfig(1) == 0
    error('Need err=[1 p] or [2 p] for error bar calculation'); 
end
if nargout == 4 && errConfig(1) == 1
    error('Cerr contains Jackknife errors: only computed when err(1) is 2'); 
end

%% Preparation
[nFreqs, nTapers, nChannels] = size(j1);
errType = errConfig(1);
pValue = errConfig(2);
probPercent = 1 - pValue/2;

%% Do the job
% --- 1. Find the number of degrees of freedom ---
if trialAve
   dim = nTapers * nChannels;
   dof = 2 * dim;
   dof1 = dof;
   dof2 = dof;
   nChannels = 1; % Collapsed due to averaging
   
   if nargin >= 6 && ~isempty(numSpikes1) 
      totSpikes1 = sum(numSpikes1);
      dof1 = fix(2 * totSpikes1 * dof / (2 * totSpikes1 + dof));
   end
   if nargin == 7 && ~isempty(numSpikes2)
      totSpikes2 = sum(numSpikes2);
      dof2 = fix(2 * totSpikes2 * dof / (2 * totSpikes2 + dof));
   end
   dof = min(dof1, dof2);
   
   % Reshape J to combine tapers and channels for processing
   j1 = reshape(j1, nFreqs, dim);
   j2 = reshape(j2, nFreqs, dim);
else
   dim = nTapers;
   dof = 2 * dim;
   dof1 = dof;
   dof2 = dof;
   
   % Handle finite size corrections per channel if needed
   for iCh = 1:nChannels
      if nargin >= 6 && ~isempty(numSpikes1)
         totSpikes1 = numSpikes1(iCh); 
         dof1 = fix(2 * totSpikes1 * dof / (2 * totSpikes1 + dof));
      end
      if nargin == 7 && ~isempty(numSpikes2)
         totSpikes2 = numSpikes2(iCh);
         dof2 = fix(2 * totSpikes2 * dof / (2 * totSpikes2 + dof));
      end
      dof(iCh) = min(dof1, dof2);
   end
end

% --- 2. Calculate Theoretical, Asymptotic Confidence Level ---
if dof <= 2
   confC = 1;
else     
   df = 1 ./ ((dof/2) - 1);
   confC = sqrt(1 - pValue.^df);
end

% --- 3. Phase Standard Deviation & Jackknife CI ---
if errType == 1
   % Theoretical Error
   totNum = nFreqs * nChannels;
   phiStd = zeros(totNum, 1); 
   coherenceCol = reshape(coherenceMag, [totNum, 1]); 
   
   % Identify valid coherence points
   idxValid = find(abs(coherenceCol - 1) >= 1.e-16);
   
   dofMat = repmat(dof, [nFreqs, 1]);
   dofMat = reshape(dofMat, [totNum 1]); 
   
   % Compute phase error
   phiStd(idxValid) = sqrt((2 ./ dofMat(idxValid) .* (1 ./ (coherenceCol(idxValid).^2) - 1))); 
   phiStd = reshape(phiStd, [nFreqs nChannels]);
   
elseif errType == 2
    % Jackknife Error
    
    % Inverse T-distribution
    tCrit = tinv(probPercent, dof - 1);
    
    % Initialize 
    % atanhCxyk stores transformed coherence dropping 1 estimate
    % phasefactorxyk stores phase vectors
    
    for k = 1:dim % dim is the number of 'independent' estimates (tapers or tapers*trials)
        % Create indices excluding k
        idxK = setdiff(1:dim, k);
        
        j1K = j1(:, idxK, :);
        j2K = j2(:, idxK, :);
        
        % Compute spectral estimates excluding k
        eJ1K = squeeze(sum(j1K .* conj(j1K), 2));
        eJ2K = squeeze(sum(j2K .* conj(j2K), 2));
        eJ12K = squeeze(sum(conj(j1K) .* j2K, 2)); 
        
        % Compute Coherence for this drop-1 set
        cXyK = eJ12K ./ sqrt(eJ1K .* eJ2K);
        absCXyK = abs(cXyK);
        
        % Fisher Z-transform: z = sqrt(2*dim-2) * atanh(C)
        atanhCXyK(k, :, :) = sqrt(2 * dim - 2) * atanh(absCXyK); 
        phaseFactorXyK(k, :, :) = cXyK ./ absCXyK;
    end
    
    % Transform original coherence
    atanhC = sqrt(2 * dim - 2) * atanh(coherenceMag); 
    
    % Jackknife estimate of standard deviation
    % sigma = sqrt(N-1) * std(drop-1-estimates)
    sigma12 = sqrt(dim - 1) * squeeze(std(atanhCXyK, 1, 1)); 

    if nChannels == 1; sigma12 = sigma12'; end
    
    % Upper and Lower bounds in transformed space
    cUpper = atanhC + tCrit(ones(nFreqs, 1), :) .* sigma12;
    cLower = atanhC - tCrit(ones(nFreqs, 1), :) .* sigma12;
    
    % Convert back using tanh
    coherenceErr(1, :, :) = max(tanh(cLower / sqrt(2 * dim - 2)), 0); % Ensure positive
    coherenceErr(2, :, :) = tanh(cUpper / sqrt(2 * dim - 2));
    
    % Phase standard deviation via Jackknife
    phiStd = sqrt( (2 * dim - 2) * (1 - abs(squeeze(mean(phaseFactorXyK)))) );
    
    if trialAve; phiStd = phiStd'; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%