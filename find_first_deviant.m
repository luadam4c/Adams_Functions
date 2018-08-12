function [valDeviant, idxDeviant] = find_first_deviant (vector, varargin)
%% Finds the index of the first deviant from preceding peers in a time series
% Usage: [valDeviant, idxDeviant] = find_first_deviant (vector, varargin)
% Outputs:
%       valDeviant  - value of the deviant 
%                       (large difference between consecutive points)
%                   specified as a numeric scalar
%       idxDeviant  - starting index of the deviant
%                   specified as a positive integer scalar
%
% Arguments:    
%       vector      - a time series
%                   must be a numeric vector
%       varargin    - 'Deviant2Peers': deviant to peers ratio
%                   must be a positive scalar
%                   default == 2
%                   - 'PeersWindowSize': peers window size in samples
%                   must be a positive integer scalar
%                   default == 5 samples
%
%
% Requires:
%       
% Used by:
%       /home/Matlab/Adams_Functions/find_first_jump.m

% File History:
% 2018-08-11 Adapted from find_first_jump.m

%% Default values for optional arguments
deviant2peersDefault = 2;        % default deviant to peers ratio
peersWindowSizeDefault = 5;     % default peers window in samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vector', ...                  % voltage vector
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Deviant2peers', deviant2peersDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'PeersWindowSize', peersWindowSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, vector, varargin{:});
deviant2peers = iP.Results.Deviant2peers;
peersWindowSize = iP.Results.PeersWindowSize;

%% Do the job
% Get the number of deviants
nDeviants = length(vector);

% Get the starting index of the first deviant to test
idxFirstDeviantToTest = peersWindowSize + 1;

% Iterate through all peers windows
for idxDeviant = idxFirstDeviantToTest:nDeviants
    % Get the indices of the peers window
    indPeers = (idxDeviant - peersWindowSize):(idxDeviant - 1);

    % Compute the root-mean-square average of the peers
    avgPeers = rms(vector(indPeers));

    % Get the value of the deviant
    valDeviant = vector(idxDeviant);

    % If the deviant magnitude is more than deviant2peers of avgPeers, stop
    if abs(valDeviant) > avgPeers
        break
    elseif idxDeviant == nDeviants
        idxDeviant = [];
        valDeviant = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:


%}
