function [output1] = plot_evoked_LFP (fileName, varargin)
%% Plots an evoked local field potential
% Usage: [output1] = plot_evoked_LFP (fileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:    
%       fileName     - TODO: Description of fileName
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/parse_abf.m
%
% Used by:    
%       /TODO:dir/TODO:file

% File History:
% 2018-09-17 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'fileName', ...                  % TODO: Description of fileName
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, fileName, varargin{:});
param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Load data and prepare for plotting
% Load and parse the abf file
[data, siUs, abfParams, tVec, vVecs, iVecs] = parse_abf(fileName);

% Extract the parsed parameters
channelTypes = abfParams.channelTypes;
channelUnits = abfParams.channelUnits;
channelLabels = abfParams.channelLabels;
nDimensions = abfParams.nDimensions;
nSamples = abfParams.nSamples;
nChannels = abfParams.nChannels;
nSweeps = abfParams.nSweeps;

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
