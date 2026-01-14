function [jData, meanSpikes, nSpikes] = mtfftpt (data, tapers, nFft, timeGrid, freqs, freqIndices)
%% Multi-taper fourier transform for point process given as times
% Usage: [jData, meanSpikes, nSpikes] = mtfftpt (data, tapers, nFft, timeGrid, freqs, freqIndices)
% Explanation:
%       Computes the Fourier transform of a point process (e.g., spike times)
%       using multi-taper estimation.
%
% Outputs:
%       jData       - FFT 
%                   (form: frequency index x taper index x channels/trials)
%       meanSpikes  - Number of spikes per sample in each channel
%       nSpikes     - Total number of spikes in each channel
%
% Arguments:
%       data        - Point process data
%                   (struct array of times with dimension channels/trials; 
%                    also takes in 1d array of spike times as a column vector) 
%       tapers      - Precalculated tapers from dpss
%       nFft        - Length of padded data
%       timeGrid    - Time points at which tapers are calculated
%       freqs       - Frequencies of evaluation
%       freqIndices - Indices corresponding to frequencies f
%
% Requires:
%
% Used by:
%       coherencycpt.m
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from mtfftpt.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 6
    error('Need all input arguments'); 
end

%% Preparation
if isstruct(data)
    nChannels = length(data); 
else
    nChannels = 1; 
end

nTapers = size(tapers, 2); 
nFreq = length(freqs); 

if nFreq ~= length(freqIndices)
    error('frequency information (last two arguments) inconsistent'); 
end

% FFT of tapers
hTapers = fft(tapers, nFft, 1);  
hTapers = hTapers(freqIndices, :); % restrict to required frequencies

w = 2 * pi * freqs; % angular frequencies

nSpikes = zeros(1, nChannels); 
meanSpikes = zeros(1, nChannels);

%% Do the job
for ch = 1:nChannels
    % Extract data for current channel
    if isstruct(data)
        fNames = fieldnames(data);
        % Use dynamic field reference
        dataTemp = data(ch).(fNames{1});
        
        % Filter spikes within timeGrid range
        indices = find(dataTemp >= min(timeGrid) & dataTemp <= max(timeGrid));
        if ~isempty(indices)
            dataTemp = dataTemp(indices);
        end
    else
        dataTemp = data;
        indices = find(dataTemp >= min(timeGrid) & dataTemp <= max(timeGrid));
        if ~isempty(indices)
            dataTemp = dataTemp(indices);
        end
    end
    
    nSpikes(ch) = length(dataTemp);
    meanSpikes(ch) = nSpikes(ch) / length(timeGrid);
    
    if meanSpikes(ch) ~= 0
        % Interpolate tapers to spike times
        dataProj = interp1(timeGrid', tapers, dataTemp);
        
        % Compute exponential component
        exponential = exp(-1i * w' * (dataTemp - timeGrid(1))');
        
        % Compute FFT component for this channel
        jData(:, :, ch) = exponential * dataProj - hTapers * meanSpikes(ch);
    else
        jData(1:nFreq, 1:nTapers, ch) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%