function [idxPeaks, valPeaks] = adjust_peaks (data, idxPeaksAppr, directionFactor, adjustWindow, varargin)
%% Adjusts peak indices and values given approximate peak indices
% Usage: [idxPeaks, valPeaks] = adjust_peaks (data, idxPeaksAppr, directionFactor, adjustWindow, varargin)
% Explanation:
%   Given a direction factor (1 = finding maximum and -1 = finding minimum)
%   approximate peak indices and a window to look on each side, determine 
%   the actual peak indices and values
% Outputs:
%       idxPeaks    - actual peak indices
%                   specified as a positive integer vector
%       valPeaks    - actual peak values
%                   specified as a numeric vector
% Arguments:    
%       data            - raw data
%                       must be a numeric vector
%       idxPeaksAppr    - approximate peak indices
%                       must be a positive integer vector
%       directionFactor - direction factor:
%                           1 = finding maximum
%                           -1 = finding minimum
%                       must be a 1 or -1
%       adjustWindow    - adjustment window (in samples)  either side 
%                           to look for actual peaks
%                       must be a positive integer scalar
%       varargin        - 'LeftBounds': left boundaries (samples) for 
%                                       adjustment window
%                       must be a numeric vector
%                       default == []
%                       - 'RightBounds': right boundaries (samples) for 
%                                       adjustment window
%                       must be a numeric vector
%                       default == []
%
% Used by:    
%       /home/Matlab/Kojis_Functions/find_directional_events.m
%       cd/minEASE_gui_examine_events.m
%
% File History:
% 2017-06-05 Adapted from getRealMinis.m by Mark P Beenhakker
% 

leftBoundsDefault = [];             % default left boundaries (samples) 
                                    %   for adjustment window
rightBoundsDefault = [];            % default right boundaries (samples)
                                    %   for adjustment window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'adjust_peaks';

% Add required inputs to an Input Parser
% TODO

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LeftBounds', leftBoundsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'RightBounds', rightBoundsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
% parse(iP, data, idxPeaksAppr, directionFactor, adjustWindow, varargin{:}); TODO
parse(iP, varargin{:});
leftBounds = iP.Results.LeftBounds;
rightBounds = iP.Results.RightBounds;

%% Extract from arguments
nSamples = length(data);                % number of samples in data vector
nPeaks = length(idxPeaksAppr);          % number of peaks to adjust

%% Initialize actual peak indices and values
idxPeaks = zeros(size(idxPeaksAppr));
valPeaks = zeros(size(idxPeaksAppr));

%% Find actual peak indices and values
for iPeak = 1:nPeaks
    % Find the left endpoint of data range to examine
    left = max(idxPeaksAppr(iPeak) - adjustWindow, 1);
    if ~isempty(leftBounds) && ~isempty(leftBounds(iPeak))
        left = max(left, ceil(leftBounds(iPeak)));
    end

    % Find the right endpoint of data range to examine
    right = min(idxPeaksAppr(iPeak) + adjustWindow, nSamples);
    if ~isempty(rightBounds) && ~isempty(rightBounds(iPeak))
        right = min(right, floor(rightBounds(iPeak)));
    end

    % Find the data range to examine
    dataRange = data(left:right);

    % Find the peak index in the data range
    [~, idxPeakDataRange] = max(dataRange * directionFactor);

    % Find the peak index in the original data vector
    idxPeaks(iPeak) = left - 1 + idxPeakDataRange;

    % Find the peak index in the original data vector
    valPeaks(iPeak) = data(idxPeaks(iPeak));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
