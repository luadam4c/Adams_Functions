function surfaceArea = compute_surface_area (lengths, diameters, varargin)
%% Computes the surface area of a cylindrical compartmental model cell based on lengths and diameters
% Usage: surfaceArea = compute_surface_area (lengths, diameters, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       surfaceArea - surface area excluding ends (cm^2)
%                   specified as a nonnegative vector
% Arguments:
%       lengths     - lengths (um) for each compartment (each column)
%                   Note: Each row is computed independently
%                   must be a nonnegative array
%       diameters   - diameters (um) for each compartment (each column)
%                   Note: Each row is computed independently
%                   must be a nonnegative array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       /TODO:dir/TODO:file
%
% Used by:
%       cd/m3ha_create_initial_neuronparams.m

% File History:
% 2018-11-13 Adapted from ~/m3ha/optimizer4gabab/singleneuronfitting41.m
% 

%% Hard-coded constants
CM_PER_UM = 1e-4;

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
addRequired(iP, 'lengths', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'diameters', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, lengths, diameters, varargin{:});

%% Do the job
% Compute the area of the model cell in cm^2
surfaceArea = sum(pi .* diameters .* lengths .* (CM_PER_UM)^2, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%