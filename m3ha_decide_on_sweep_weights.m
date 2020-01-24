function sweepWeights = m3ha_decide_on_sweep_weights (sweepWeights, fileNames, varargin)
%% Set default weights for fitting for the GAT blockade project
% Usage: sweepWeights = m3ha_decide_on_sweep_weights (sweepWeights, fileNames, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       sweepWeights    - TODO: Description of sweepWeights
%                       must be a TODO
%       fileNames       - TODO: Description of fileNames
%                       must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_rank_neurons.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting99.m and beyond

% File History:
% 2020-01-24 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'sweepWeights');
addRequired(iP, 'fileNames');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, sweepWeights, fileNames, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if isempty(sweepWeights)
    % Count the number of files
    nFiles = numel(fileNames);

    % Weight the later files more if nFiles is a multiple of 4
    if mod(nFiles, 4) ~= 0
        sweepWeights = ones(nFiles, 1);
    else
        % Count the number of gIncr conditions to fit
        nGIncr = nFiles / 4;

        % Make sweeps with higher gIncrs weighted more
        sweepWeights = repmat(transpose(1:nGIncr), 4, 1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%