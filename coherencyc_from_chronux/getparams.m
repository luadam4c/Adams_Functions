function [tapers, pad, fs, fPass, err, trialAve, params] = getparams (params)
%% Helper function to convert structure params to variables used by Chronux routines
% Usage: [tapers, pad, fs, fPass, err, trialAve, params] = getparams (params)
% Explanation:
%       Checks the input 'params' structure for required fields.
%       Fills in default values if fields are missing or empty.
%       Converts taper parameters from [W T p] form to [TW K] form if necessary.
%
% Outputs:
%       tapers      - Tapers in [TW K] form or precalculated tapers
%       pad         - Padding factor for the FFT
%       fs          - Sampling frequency
%       fPass       - Frequency band to be used
%       err         - Error calculation mode
%       trialAve    - Flag for trial averaging
%       params      - The updated structure
%
% Arguments:
%       params      - Structure with fields:
%           .tapers : precalculated tapers from dpss or in one of the following forms:
%                       (1) A numeric vector [TW K] where TW is the
%                           time-bandwidth product and K is the number of
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
%                       Default: [3 5] (TW=3, K=5)
%           .pad    : (padding factor for the FFT) - optional
%                       Can take values -1, 0, 1, 2...
%                       -1 corresponds to no padding.
%                       0 corresponds to padding to the next highest power of 2 etc.
%                       e.g. For N = 500:
%                           if PAD = -1, we do not pad;
%                           if PAD = 0, we pad the FFT to 512 points;
%                           if pad = 1, we pad to 1024 points etc.
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
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from getparams.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Do the job

% --- Tapers ---
if ~isfield(params, 'tapers') || isempty(params.tapers)
     disp('tapers unspecified, defaulting to params.tapers=[3 5]');
     params.tapers = [3 5];
end

% Check if tapers are in form [W T p] (length 3) and convert to [TW K]
if ~isempty(params) && length(params.tapers) == 3 
    % Compute time-bandwidth product (TW = T * W)
    % params.tapers(1) is W (bandwidth)
    % params.tapers(2) is T (duration)
    timeBandwidth = params.tapers(2) * params.tapers(1);
    
    % Compute number of tapers (K = 2TW - p)
    % params.tapers(3) is p
    numTapers = floor(2 * timeBandwidth - params.tapers(3));
    
    params.tapers = [timeBandwidth numTapers];
end

% --- Padding ---
if ~isfield(params, 'pad') || isempty(params.pad)
    params.pad = 0;
end

% --- Sampling Frequency ---
if ~isfield(params, 'Fs') || isempty(params.Fs)
    params.Fs = 1;
end

% --- Frequency Band ---
if ~isfield(params, 'fpass') || isempty(params.fpass)
    params.fpass = [0 params.Fs/2];
end

% --- Error Calculation ---
if ~isfield(params, 'err') || isempty(params.err)
    params.err = 0;
end

% --- Trial Averaging ---
if ~isfield(params, 'trialave') || isempty(params.trialave)
    params.trialave = 0;
end

%% Output results
tapers   = params.tapers;
pad      = params.pad;
fs       = params.Fs;
fPass    = params.fpass;
err      = params.err;
trialAve = params.trialave;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%