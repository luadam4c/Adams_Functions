function specErr = specerr (spectrum, jData, errConfig, trialAve, numSpikes)
%% Function to compute lower and upper confidence intervals on the spectrum 
% Usage: specErr = specerr (spectrum, jData, errConfig, trialAve, numSpikes)
% Explanation:
%       Computes error bars for spectra using either asymptotic (theoretical)
%       estimates or Jackknife estimates.
%
% Outputs:
%       specErr     - Error estimates. 
%                   specErr(1, ...) is lower confidence level.
%                   specErr(2, ...) is upper confidence level.
%
% Arguments:
%       spectrum    - Spectrum data
%       jData       - Tapered fourier transforms 
%       errConfig   - [errType pValue] 
%                   errType=1: asymptotic estimates
%                   errType=2: Jackknife estimates
%                   pValue: p value for error estimates
%       trialAve    - 0: no averaging; 1: perform trial averaging
%       numSpikes   - (Optional) Number of spikes in each channel. 
%                   Specify only when finite size correction is required 
%                   (point process data only).
%
% Requires:
%       chi2inv (MATLAB Statistics Toolbox)
%       tinv (MATLAB Statistics Toolbox)
%
% Used by:
%       mtspectrumc.m
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from specerr.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 4
    error('Need at least 4 input arguments'); 
end
if errConfig(1) == 0
    error('Need err=[1 p] or [2 p] for error bar calculation. Make sure you are not asking for the output of Serr'); 
end

%% Preparation
[nFreqs, nTapers, nChannels] = size(jData);
errType = errConfig(1);
pValue  = errConfig(2);
probP   = 1 - pValue/2;
probQ   = 1 - probP;

%% Do the job
% Calculate Degrees of Freedom (dof)
if trialAve
    dim = nTapers * nChannels;
    nChannels = 1; % Effectively one channel after averaging
    degOfFreedom = 2 * dim;
    if nargin == 5
        % Finite size correction for trial averaged data
        degOfFreedom = fix(1 / (1/degOfFreedom + 1/(2 * sum(numSpikes)))); 
    end
    % Reshape J for calculation
    jData = reshape(jData, nFreqs, dim);
else
    dim = nTapers;
    degOfFreedom = 2 * dim * ones(1, nChannels);
    for ch = 1:nChannels
        if nargin == 5
            % Finite size correction per channel
            degOfFreedom(ch) = fix(1 / (1/degOfFreedom(ch) + 1/(2 * numSpikes(ch)))); 
        end 
    end
end

specErr = zeros(2, nFreqs, nChannels);

% Calculate Error Bars based on Type
if errType == 1
    % --- Asymptotic Estimates ---
    invChiSqP = chi2inv(probP, degOfFreedom);
    invChiSqQ = chi2inv(probQ, degOfFreedom);
    
    % Expand matrices for element-wise operation
    specErr(1, :, :) = degOfFreedom(ones(nFreqs, 1), :) .* spectrum ./ invChiSqP(ones(nFreqs, 1), :);
    specErr(2, :, :) = degOfFreedom(ones(nFreqs, 1), :) .* spectrum ./ invChiSqQ(ones(nFreqs, 1), :);
    
elseif errType == 2
    % --- Jackknife Estimates ---
    tCrit = tinv(probP, dim - 1);
    
    % 1-drop projection loop
    for k = 1:dim
        indices = setdiff(1:dim, k);
        jJackknife = jData(:, indices, :); 
        eJJackknife = squeeze(sum(jJackknife .* conj(jJackknife), 2));
        sJackknife(k, :, :) = eJJackknife / (dim - 1); % 1-drop spectrum
    end
    
    % Calculate standard deviation of log spectrum
    sigma = sqrt(dim - 1) * squeeze(std(log(sJackknife), 1, 1)); 
    
    if nChannels == 1; sigma = sigma'; end
    
    conf = repmat(tCrit, nFreqs, nChannels) .* sigma;
    conf = squeeze(conf); 
    
    specErr(1, :, :) = spectrum .* exp(-conf); 
    specErr(2, :, :) = spectrum .* exp(conf);
end

%% Output results
specErr = squeeze(specErr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%