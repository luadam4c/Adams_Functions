function gpas = compute_gpas (inputResistance, surfaceArea, varargin)
%% Computes the passive conductance (gpas, in S/cm^2) from input resistance and surface area
% Usage: gpas = compute_gpas (inputResistance, surfaceArea, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       gpas        - passive conductance (S/cm^2)
%                   specified as a nonnegative vector
% Arguments:
%       inputResistance - input resistance (MOhm)
%                       must be a nonnegative vector
%       surfaceArea     - surface area (cm^2)
%                       must be a nonnegative vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:
%       cd/m3ha_create_initial_neuronparams.m

% File History:
% 2018-11-13 Adapted from ~/m3ha/optimizer4gabab/singleneuronfitting41.m
% 

%% Hard-coded constants
OHM_PER_MOHM = 1e6;

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
addRequired(iP, 'inputResistance', ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'vector'}));
addRequired(iP, 'surfaceArea', ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'vector'}));

% Read from the Input Parser
parse(iP, inputResistance, surfaceArea, varargin{:});

%% Do the job
% [S/cm^2] = 1 / ([MOhm] * [Ohm / MOhm] * [cm^2])
gpas = 1 ./ (inputResistance .* OHM_PER_MOHM .* surfaceArea);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%