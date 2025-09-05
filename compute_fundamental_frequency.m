function fundamentalFrequencies = compute_fundamental_frequency (vectors, varargin)
%% Computes the fundamental frequency for one or more time-series vectors
% Usage: fundamentalFrequencies = compute_fundamental_frequency (vectors, varargin)
% Explanation:
%       This function calculates the fundamental frequency of a signal using 
%       the Fast Fourier Transform (FFT). It identifies the frequency component 
%       with the highest power, excluding the DC component (0 Hz). The function
%       can process a single vector, a matrix of vectors (column-wise), or a
%       cell array of vectors.
%
% Example(s):
%       % --- Example 1: A simple 10 Hz sine wave ---
%       siMs = 1;                           % 1 ms sampling interval (1000 Hz)
%       fs = 1000 / siMs;                   % Sampling frequency in Hz
%       t = (0:siMs:1000-siMs)' ./ 1000;    % Time vector of 1 second
%       y = sin(2 * pi * 10 * t);           % 10 Hz sine wave
%       freq = compute_fundamental_frequency(y, 'SamplingIntervalMs', siMs)
%
%       % --- Example 2: Using a time vector instead ---
%       freq_from_timevec = compute_fundamental_frequency(y, 'TimeVecs', t .* 1000)
%
%       % --- Example 3: A cell array with two signals (5 Hz and 20 Hz) ---
%       y2 = cos(2 * pi * 20 * t);
%       signals = {y, y2};
%       freqs = compute_fundamental_frequency(signals, 'SamplingIntervalMs', siMs)
%       
% Outputs:
%       fundamentalFrequencies - The detected fundamental frequency in Hz for
%                                each input vector.
%                                specified as a column vector
%
% Arguments:    
%       vectors     - Time-series data vector(s).
%                   Note: If a cell array, each element must be a vector.
%                         If a numeric array, each column is treated as a vector.
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'SamplingIntervalMs': sampling interval(s) in milliseconds.
%                   Note: Only one of 'SamplingIntervalMs' or 'TimeVecs'
%                         can be provided.
%                   must be a positive scalar or a vector with one element per input vector
%                   default == 0.1
%                   - 'TimeVecs': corresponding time vectors in milliseconds.
%                   Note: If provided, the sampling interval is calculated as the
%                         average difference between consecutive time points.
%                         Only one of 'SamplingIntervalMs' or 'TimeVecs' can be provided.
%                   must be a numeric array or a cell array of numeric vectors
%                   default == []
%
% Requires:
%       cd/array_fun.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_format_vector_sets.m
%
% Used by:    
%       cd/parse_oscillation.m

% File History:
% 2025-09-04 Created by Gemini, based on logic from virt_moore.m
% 2025-09-04 Modified by Gemini to use array_fun instead of a for loop
% 

%% Hard-coded parameters
% None

%% Default values for optional arguments
samplingIntervalMsDefault = 0.1;    % default sampling interval is 0.1 ms (10 kHz)
timeVecsDefault = [];               % no time vectors by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;                                                   % Create an input parser object
iP.FunctionName = mfilename;                                        % Set the function name for error messages

% Add required inputs to the Input Parser
addRequired(iP, 'vectors', ...                                      % The input signal(s)
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...               % Validation function to ensure input is numeric
                'vectors must be a numeric array or a cell array of numeric vectors!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SamplingIntervalMs', samplingIntervalMsDefault, ... % The sampling interval in ms
    @(x) validateattributes(x, {'numeric'}, {'positive'}));         % Must be a positive number

addParameter(iP, 'TimeVecs', timeVecsDefault, ...                   % The time vector(s) in ms
    @(x) isempty(x) || isnumeric(x) || iscellnumeric(x));           % Can be empty, numeric, or a cell array

% Read from the Input Parser
parse(iP, vectors, varargin{:});                                    % Parse the inputs
siMs = iP.Results.SamplingIntervalMs;                               % Store the sampling interval
timeVecs = iP.Results.TimeVecs;                                     % Store the time vectors

% Check that 'SamplingIntervalMs' and 'TimeVecs' are not both provided
if ~isempty(timeVecs) && ~any(siMs ~= samplingIntervalMsDefault)
    % If TimeVecs is provided, clear the default siMs to prioritize TimeVecs
    siMs = [];                                                      
elseif ~isempty(timeVecs) && any(siMs ~= samplingIntervalMsDefault)
    % If both are explicitly provided, raise an error
    error('Usage error: ''SamplingIntervalMs'' and ''TimeVecs'' cannot both be provided.');
end

%% Preparation
% Force vectors to be a column cell array for consistent processing
vectors = force_column_cell(vectors);                               

% Determine the number of vectors to process
nVectors = numel(vectors);                                          

% Handle time information
if ~isempty(timeVecs)
    % If time vectors are provided, force them into a column cell array format
    timeVecs = force_column_cell(timeVecs);                         

    % Ensure the number of time vectors matches the number of data vectors
    [timeVecs, vectors] = match_format_vector_sets(timeVecs, vectors);

    % Calculate the sampling interval for each vector from its time vector
    % The sampling interval is the mean of the differences between consecutive time points
    siMs = cellfun(@(v) mean(diff(v)), timeVecs);                    
else
    % If a single sampling interval is given, replicate it for all vectors
    if isscalar(siMs)                                               
        siMs = repmat(siMs, nVectors, 1);                           
    end
end

%% Do the job
% Compute fundamental frequency for all vectors using array_fun
% This approach passes each vector and its corresponding sampling interval
% to the local function 'compute_single_freq' for processing.
fundamentalFrequencies = array_fun(@compute_single_freq, vectors, num2cell(siMs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freqFundamental = compute_single_freq(thisVector, thisSiMs)
%% Computes the fundamental frequency for a single vector

% Return NaN if the vector is empty or too short for FFT
if isempty(thisVector) || numel(thisVector) < 2
    freqFundamental = NaN;
    return;
end

% --- Begin FFT Analysis (adapted from virt_moore.m) ---

% Calculate sampling frequency (fs) in Hz from sampling interval in ms
fs = 1000 / thisSiMs;

% Get the length of the signal (number of samples)
L = length(thisVector);

% Perform the Fast Fourier Transform on the signal
Y = fft(thisVector);

% Compute the two-sided power spectrum P2
P2 = abs(Y/L);

% Compute the single-sided spectrum P1 from the first half of P2
P1 = P2(1:floor(L/2)+1);

% Double the amplitudes in P1 to account for the energy from the negative frequencies
% Skip the first (DC) and last (Nyquist) components
P1(2:end-1) = 2*P1(2:end-1);

% Create the frequency vector corresponding to the P1 spectrum
fVec = fs*(0:(L/2))/L;

% Find the index of the maximum power in the spectrum, EXCLUDING the DC component (P1(1))
[~, idx] = max(P1(2:end));

% Get the fundamental frequency from the frequency vector
% We add 1 to the index because the search started from the 2nd element of P1
freqFundamental = fVec(idx + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

