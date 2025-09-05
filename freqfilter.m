function filteredData = freqfilter (data, fc, varargin)
%% Uses a Butterworth filter twice to filter data by a frequency band (each column is a vector of samples)
% Usage: filteredData = freqfilter (data, fc, si (opt), varargin)
% Explanation: 
%       A Butterworth filter is applied twice (once forward, once backward)
%       to remove lag effect of filter.
%
% Outputs:
%       filteredData - The filtered version of data. The output format (numeric
%                      matrix or cell array) matches the input format of 'data'.
%
% Arguments:    
%       data        - Data to be filtered. Each column is a vector of samples.
%                   Note: If 'data' is a cell array of numeric vectors, the
%                         output will also be a cell array.
%                   must be a nonempty numeric array or a cell array of numeric vectors
%       fc          the cutoff frequency(ies) (Hz or normalized) for the filter
%                   must be a numeric and:
%                       a scalar by default or if ftype == 'low' or 'high'
%                       a two-element vector if ftype == 'bandpass' or 'stop'
%                   consistent with the documentation for butter()
%       si          - (opt) sampling interval (seconds)
%                   must be a numeric scalar
%                   default == 1 (fc must have normalized units in this case)
%       varargin    - 'FilterType': filter type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'low'       - lowpass filter with cutoff frequency fc
%                       'high'      - highpass filter with cutoff frequency fc
%                       'bandpass'  - bandpass filter of order 2*npoles 
%                                       with cutoff frequencies fc(1) & fc(2)
%                       'stop'      - bandstop filter of order 2*npoles 
%                                       with cutoff frequencies fc(1) & fc(2)
%                   default == 'low' if fc has one element and 
%                           == 'bandpass' if fc has two elements
%                   - 'FilterOrder': order of filter
%                       i.e., the number of poles in the transfer function
%                       i.e., the order of the polynomial in the denominator 
%                           of the transfer function
%                       Note: the higher the order the steeper the cutoff
%                   must be a numeric scalar that is a positive integer
%                   consistent with the documentation for butter()
%
% Requires:
%
% Used by:
%       cd/crosscorr_profile.m
%       cd/detect_spikes_multiunit.m
%       cd/filter_and_extract_pulse_response.m
%       cd/identify_CI_protocol.m
%       cd/parse_oscillation.m
%       cd/parse_pleth_trace.m
%       cd/plot_calcium_imaging_traces.m
%       cd/minEASE.m
%           many others; apply this command in a LINUX terminal to find them:
%             grep --include=*.m -rlw '/home/Matlab/' -e "freqfilter"
%       \Shared\Code\vIRt\virt_moore.m
%
% File History:
% 2018-08-03 AL - Adapted from /home/Matlab/Marks_Functions/zof_mark.m
% 2018-08-03 Made npoles an optional parameter 'FilterOrder'
% 2018-08-03 Made si an optional parameter
% 2020-08-13 Improved bandpass filter design based on 
%               https://dsp.stackexchange.com/questions/42797/why-butterworth-bandpass-filter-changes-signal-in-the-passband
% 2021-05-07 Added maxNPoles and set it at 20
% 2025-09-04 Improved by Gemini to handle cell arrays and check fc bounds.

%% Hard-coded parameters
validFilterTypes = {'low', 'high', 'bandpass', 'stop', 'auto'};
passBandRipple = 3;
stopBandAttenuation = 40;
passStopBandDiffHz = 0.5;
maxNPoles = 20;

%% Default values for optional arguments
defaultSamplingInterval = [];       % to be set later
defaultFilterType = 'auto';         % to be set later
defaultFilterOrder = [];            % use an 8th order Butterworth filter
                                    %   by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'data', ...     % vector of samples
    @(x) assert(isnumeric(x) || iscell(x), ...
                'Data must be a numeric array or a cell array of numeric vectors!'));
addRequired(iP, 'fc', ...       % the cutoff frequency(ies) (Hz or normalized)
    @(x) isnumeric(x) && isvector(x) && numel(x) <= 2);

% Add optional inputs to the Input Parser
addOptional(iP, 'si', defaultSamplingInterval, ... % sampling interval (seconds)
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FilterType', defaultFilterType, ...       % filter type
    @(x) any(validatestring(x, validFilterTypes)));
addParameter(iP, 'FilterOrder', defaultFilterOrder, ...     % order of filter
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, data, fc, varargin{:});
si = iP.Results.si;
ftype = validatestring(iP.Results.FilterType, validFilterTypes);
nPoles = iP.Results.FilterOrder;

% Set dependent argument defaults
if isempty(si)
    % In this case fc must have normalized units
    if max(fc) >= 1 || min(fc) <= 0
        error(['A sampling interval (seconds) must be provided', ...
                ' if cutoff frequencies are not normalized!']);
    end
    
    % Set si to 1
    si = 1; 
end

if strcmp(ftype, 'auto')
    if numel(fc) > 1
        ftype = 'bandpass';
    else
        ftype = 'low';
    end
end

%% Preparation
% Store whether the input was a cell array
isCellInput = iscell(data);

% Calculate the Nyquist frequency
nyquistFreq = 1 / (2 * si);

% If fc has two elements but ftype is not 'bandpass' or 'stop',
%   return error
if numel(fc) > 1 && ~any(strcmp(ftype, {'bandpass', 'stop'}))
    error(['Filter type must be ''bandpass'' or ''stop'' ', ...
            'if two cutoff frequencies are provided!']);
end

% Validate fc against lower and upper bounds
if any(fc < 0)
    warning('Cutoff frequency has negative values, will change to 0!');
    fc(fc < 0) = 0;
end
if any(fc >= nyquistFreq)
    warning('Cutoff frequency has values greater than Nyquist frequency (%.2f Hz).', nyquistFreq);
    warning('Will change to 0.99 times Nyquist frequency!');
    fc(fc >= nyquistFreq) = 0.99 * nyquistFreq;
end

% Change bandpass to lowpass if necessary
if strcmp(ftype, 'bandpass') && fc(1) <= 0
    ftype = 'low';
    fc = fc(2);
end

% Change bandstop to highpass if necessary
if strcmp(ftype, 'stop') && fc(1) <= 0
    ftype = 'high';
    fc = fc(2);
end

% Prepare data and get dimensions if input is a matrix
if ~isCellInput
    % If data is a vector, make sure it is a column
    if isvector(data)
        data = data(:);
    end

    % Count the number of traces to filter
    nTraces = size(data, 2);

    % Count the number of samples per trace
    nSamples = size(data, 1);
end

% Force fc as a column vector
fc = fc(:);

% Compute filter order and cutoff frequencies
if isempty(nPoles)
    % Compute the pass band frequencies
    passBandFrequency = fc / nyquistFreq;

    % Compute the stop band frequencies
    switch ftype
        case 'low'
            stopBandFrequency = (fc + passStopBandDiffHz) / nyquistFreq;
        case 'high'
            stopBandFrequency = (fc - passStopBandDiffHz) / nyquistFreq;
        case 'bandpass'
            stopBandFrequency = (fc + passStopBandDiffHz * [-1; 1]) / nyquistFreq;
        case 'stop'
            stopBandFrequency = (fc + passStopBandDiffHz * [1; -1]) / nyquistFreq;
    end

    % Fix values if out of range
    stopBandFrequency(stopBandFrequency >= 1) = ...
        (passBandFrequency(stopBandFrequency >= 1) + 1) / 2;
    stopBandFrequency(stopBandFrequency <= 0) = ...
        (passBandFrequency(stopBandFrequency <= 0) + 0) / 2;
    
    % Compute an appropriate filter order if not provided
    [nPoles, Wn] = buttord(passBandFrequency, stopBandFrequency, ...
                            passBandRipple, stopBandAttenuation);
                        
    % Prevent the filter order from exceeding maxNPoles
    if nPoles > maxNPoles
        nPoles = maxNPoles;
    end
else
    % Find the normalized cutoff frequency(ies) Wn = fc/(fs/2), 
    %   where fs = sampling frequency (Hz) = 1/si
    %   and fs/2 is the Nyquist frequency
    Wn = fc / nyquistFreq;      % normalized cutoff frequency (half-cycles/sample)
end

% Find the transfer function coefficients of a Butterworth filter
%   with order npoles and normalized cutoff frequency(ies) Wn
[numeratorCoeff, denominatorCoeff] = butter(nPoles, Wn, ftype);

%% Filter each trace one by one
if isCellInput
    % If input is a cell array, filter each cell individually and return a cell array
    filteredData = cellfun(@(x) filter_single_trace(x, numeratorCoeff, denominatorCoeff), ...
                           data, 'UniformOutput', false);
else
    % If input is a matrix, filter each column
    
    % Check if there are enough data points for the filter
    orderFilter = filtord(numeratorCoeff, denominatorCoeff);
    if nSamples <= 3 * orderFilter
        error(['Not enough data points (%d) to apply a ', ...
                'Butterworth filter of order %d twice!\n'], ...
                nSamples, orderFilter);
    end

    % Initialize filtered data
    filteredData = zeros(size(data));

    for i = 1:nTraces
        % Extract the ith trace to filter
        thisTrace = data(:, i);

        % Ignore NaNs from padding
        isNan = isnan(thisTrace);
        thisTrace(isNan) = 0; % Replace NaNs with 0 for filtering

        % Lowpass-filter data twice (forward & reverse directions)
        thisTraceFiltered = filtfilt(numeratorCoeff, denominatorCoeff, thisTrace);

        % Restore NaNs
        thisTraceFiltered(isNan) = NaN;

        % Place filtered trace in output matrix
        filteredData(:, i) = thisTraceFiltered;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outTrace = filter_single_trace(inTrace, b, a)
%% Applies filtfilt to a single vector, with safety checks.

% Get the effective order of the filter for the length check
orderFilter = filtord(b, a);

% If trace is empty or too short for the filter, return it unchanged
if isempty(inTrace) || numel(inTrace) <= 3 * orderFilter
    outTrace = inTrace;
    warning(['Not enough data points to apply a ', ...
            'Butterworth filter of order %d twice!'], ...
            orderFilter);
    warning('The original vector will be returned!');
    return;
end

% Apply the filter forwards and backwards
outTrace = filtfilt(b, a, inTrace);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

function filtered_traces = freqfilter(data,lp,poles,si)
% data is a vector of samples
% lp is the low pass cutoff
% poles is number of poles?
% not so clear on meaning of "poles" here
% si is sample interval in sec
% output is filtered trace that was run through a butterworth
%filter twice, once forward, then once backward to remove lag effect of
%filter
filtered_traces = zeros(size(data,1), size(data,2)) ;
filtered_traces(:,1) = data(:,1) ;

for i = 1:size(data,2)
sweep = data(:,i) ;
[blpf,alpf]=butter(poles,lp*2*si);
filttr=filter(blpf,alpf,sweep);
filttr=filter(blpf,alpf,filttr(end:-1:1)); % filter the reversed trace to remove time lag from filter
filttr=filttr(end:-1:1); % reverse the trace to original orientation
filtered_traces(:,i) = filttr ;

end

% Extract the ith trace to filter
thisTrace = data(:, i);             % ith trace to filter

% Find the normalized cutoff frequency Wn = fc/(fs/2), 
%   where fs = sampling frequency (Hz) = 1/si 
%   and fs/2 is the Nyquist frequency
Wn = fc * 2 * si;       % normalized cutoff frequency (half-cycles/sample)

% Find the transfer function coefficients of a Butterworth filter
%   with order npoles and normalized cutoff frequency Wn
[numeratorCoeff, denominatorCoeff] = butter(npoles, Wn, ftype);    

% Apply filter to trace
filtTr1 = filter(numeratorCoeff, denominatorCoeff, thisTrace);

% Filter the reversed trace to remove time lag from filter
filtTr2 = filter(numeratorCoeff, denominatorCoeff, filtTr1(end:-1:1)); 

% Reverse the trace to original orientation
filtTr3 = filtTr2(end:-1:1);

% Place filtered trace in output matrix
filteredData(:, i) = filtTr3;

@(x) isnumeric(x) && isvector(x) && numel(x) <= 2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
