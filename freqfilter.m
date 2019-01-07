function filteredData = freqfilter (data, fc, varargin)
%% Uses a Butterworth filter twice to filter data by a frequency band (each column is a vector of samples)
% Usage: filteredData = freqfilter (data, fc, varargin)
% Explanation: 
%       A Butterworth filter is applied twice (once forward, once backward)
%       to remove lag effect of filter.
% Outputs:
%       filteredData - the filtered version of data
% Arguments:    
%       data        - data where each column is a vector of samples
%                   must be a nonempty numeric array
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
%                   - 'NPoles': order of filter
%                       i.e., the number of poles in the transfer function
%                       i.e., the order of the polynomial in the denominator 
%                           of the transfer function
%                       Note: the higher the order the steeper the cutoff
%                   must be a numeric scalar that is a positive integer
%                   consistent with the documentation for butter()
%
% Used by:
%       cd/filter_and_extract_pulse_response.m
%       cd/identify_CI_protocol.m
%       /home/Matlab/minEASE/minEASE.m
%           many others; apply this command in a LINUX terminal to find them:
%             grep --include=*.m -rlw '/home/Matlab/' -e "freqfilter"
%
% File History:
% 2018-08-03 AL - Adapted from /home/Matlab/Marks_Functions/zof_mark.m
% 2018-08-03 Made npoles an optional parameter 'FilterOrder'
% 2018-08-03 Made si an optional parameter
% TODO: Check lower and upper bounds for fc (0, Nyquist frequency)
% TODO: Allow data to be a cell array but attempt to concatenate
%           into array 

%% Hard-coded parameters
validFilterTypes = {'low', 'high', 'bandpass', 'stop', 'auto'};

%% Default values for optional arguments
defaultSamplingInterval = [];       % to be set later
defaultFilterType = 'auto';         % to be set later
defaultFilterOrder = 8;             % use an 8th order Butterworth filter
                                    %   by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'data', ...     % vector of samples
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
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
% If fc has two elements but ftype is not 'bandpass' or 'stop',
%   return error
if numel(fc) > 1 && ~any(strcmp(ftype, {'bandpass', 'stop'}))
    error(['Filter type must be ''bandpass'' or ''stop'' ', ...
            'if two cutoff frequencies are provided!']);
end

% If data is a vector, make sure it is a column
if isvector(data)
    data = data(:);
end

% Count the number of traces to filter
nTraces = size(data, 2);

% Count the number of samples per trace
nSamples = size(data, 1);

% Find the normalized cutoff frequency(ies) Wn = fc/(fs/2), 
%   where fs = sampling frequency (Hz) = 1/si
%   and fs/2 is the Nyquist frequency
Wn = fc * 2 * si;           % normalized cutoff frequency (half-cycles/sample)

% Find the transfer function coefficients of a Butterworth filter
%   with order npoles and normalized cutoff frequency(ies) Wn
[numeratorCoeff, denominatorCoeff] = butter(nPoles, Wn, ftype);

% Check the order of the filter
orderFilter = filtord(numeratorCoeff, denominatorCoeff);
if nSamples <= 3 * orderFilter
    error(['Not enough data points to apply a ', ...
            'Butterworth filter of order %d twice!\n'], ...
            orderFilter);
end

%% Filter each trace one by one
filteredData = zeros(size(data));       % initialize filtered data
for i = 1:nTraces
    % Extract the ith trace to filter
    thisTrace = data(:, i);             % ith trace to filter

    % Lowpass-filter data twice (forward & reverse directions)
    thisTraceFiltered = filtfilt(numeratorCoeff, denominatorCoeff, thisTrace);

    % Place filtered trace in output matrix
    filteredData(:, i) = thisTraceFiltered;
end

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