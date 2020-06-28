function [parsedParams, parsedData] = parse_psd (dataVec, varargin)
%% Parses the power spectral density and compute peak frequencies of a data vector
% Usage: [parsedParams, parsedData] = parse_psd (dataVec, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       parsedParams    - TODO: Description of parsedParams
%                       specified as a TODO
%       parsedData      - TODO: Description of parsedData
%                       specified as a TODO
%
% Arguments:
%       dataVec     - TODO: Description of dataVec
%                   must be a TODO
%       varargin    - 'SamplingFrequencyHz': sampling frequency in Hz
%                   must be a positive scalar
%                   default == 1 Hz
%                   - 'FilterWindowHz': filter window for PSD in Hz
%                   must be a positive scalar
%                   default == 2 Hz
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/parse_atf_swd.m
%       ~/EEG_gui/EEG_gui.m

% File History:
% 2020-06-28 Moved from EEG_gui.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
samplingFrequencyHzDefault = 1;
filterWindowHzDefault = 2;          % 2 Hz window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'dataVec');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SamplingFrequencyHz', samplingFrequencyHzDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'FilterWindowHz', filterWindowHzDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, dataVec, varargin{:});
samplingFrequencyHz = iP.Results.SamplingFrequencyHz;
filterWindowHz = iP.Results.FilterWindowHz;

%% Do the job
% Subtract the data by its mean to get the signal
signal = dataVec - mean(dataVec);

% Signal will be padded with zeros to make length transformLength
transformLength = 2 ^ nextpow2(length(signal));

% Get the discrete Fourier transform of signal
dft = fft(signal, transformLength);

% Number of unique frequencies in Hz
nFreqs = transformLength/2 + 1;

% Unique frequencies in Hz
freqVec = samplingFrequencyHz/2 * linspace(0, 1, nFreqs);

% Part of the discrete Fourier transform of signal
%   that corresponds to unique frequencies
dftPart = dft(1:nFreqs);

% Power spectral density (power spectrum)
psd = 2 * dftPart .* conj(dftPart) * ...
        (1/transformLength) * (1/samplingFrequencyHz);

% Correct the end points (these don't need to be mirrored)
psd(1) = psd(1)/2;
psd(end) = psd(end)/2;

% Get the frequency spacing in Hz
freqSpacing = freqVec(2) - freqVec(1);

% Compute the number of frequencies within a median filter window
medianFilterWindowFreqs = round(filterWindowHz / freqSpacing);

% Take the median filter of the PSD
psdFiltered = medfilt1(psd, medianFilterWindowFreqs);

% Compute the span of a moving-average filter (must be odd)
span = 2 * floor(medianFilterWindowFreqs / 2) + 1;

% Take the moving-average filter of the median-filtered PSD
%   with the same window
psdSmoothed = smooth(psdFiltered, span);

% Use findpeaks to find local maximums
[ampPeaks, idxPeaks] = findpeaks(psdSmoothed);

% Find the number of peaks
nPeaks = length(idxPeaks);

% Find the top three frequencies
if nPeaks == 0
    % If there are no peaks, make all frequencies NaN
    freqSorted = [];
    idxSorted = [];
else
    % Sort the peaks by amplitude
    [~, peakNoSorted] = sort(ampPeaks, 'descend');

    % Select three largest peaks (indices and values) with 
    %   none of them within filterWindowHz of each other
    idxSelected = zeros(3, 1);
    freqSelected = zeros(3, 1);

    % Find the first peak
    idxSelected(1) = idxPeaks(peakNoSorted(1));
    freqSelected(1) = freqVec(idxSelected(1));

    % Find the second peak
    last1 = 1;
    if nPeaks > last1
        % Try to look for a peak that isn't within filterWindowHz of the first peak
        for iSorted = (last1 + 1):nPeaks
            idx = idxPeaks(peakNoSorted(iSorted));
            fNow = freqVec(idx);
            if abs(fNow - freqSelected(1)) > filterWindowHz
                last1 = iSorted;
                idxSelected(2) = idx;
                freqSelected(2) = fNow;
                break;
            end
        end
        % If no more peaks, make the second peak NaN
        if last1 == 1
            idxSelected(2) = NaN;
            freqSelected(2) = NaN;
        end
    else
        % If no more peaks, make the second peak NaN
        idxSelected(2) = NaN;
        freqSelected(2) = NaN;
    end

    % Find the third peak
    last2 = last1;
    if nPeaks > last2
        % Try to look for a peak that isn't within filterWindowHz 
        %   of the previous peaks
        for iSorted = (last2 + 1):nPeaks
            idx = idxPeaks(peakNoSorted(iSorted));
            fNow = freqVec(idx);
            if abs(fNow - freqSelected(1)) > filterWindowHz && ...
                abs(fNow - freqSelected(2)) > filterWindowHz
                last2 = iSorted;
                idxSelected(3) = idx;
                freqSelected(3) = fNow;
                break;
            end
        end
        % If no more peaks, make the third peak NaN
        if last2 == last1
            idxSelected(3) = NaN;
            freqSelected(3) = NaN;
        end
    else
        % If no more peaks, make the third peak NaN
        idxSelected(3) = NaN;
        freqSelected(3) = NaN;
    end

    % Order the three peaks based on frequency
    [freqSorted, origInd] = sort(freqSelected, 'ascend');
    idxSorted = idxSelected(origInd);
end

%% Output results
parsedParams.samplingFrequencyHz = samplingFrequencyHz;
parsedParams.filterWindowHz = filterWindowHz;
parsedParams.transformLength = transformLength;
parsedParams.freqSpacing = freqSpacing;

parsedData.dft = dft;
parsedData.psd = psd;
parsedData.psdFiltered = psdFiltered;
parsedData.psdSmoothed = psdSmoothed;
parsedData.freqVec = freqVec;
parsedData.freqSelected = freqSelected;
parsedData.idxSelected = idxSelected;
parsedData.freqSorted = freqSorted;
parsedData.idxSorted = idxSorted;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%