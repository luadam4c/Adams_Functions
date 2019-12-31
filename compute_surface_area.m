function surfaceArea = compute_surface_area (lengths, diameters, varargin)
%% Computes the surface area(s) (cm^2) of cylindrical compartmental model cell(s) based on lengths (um) and diameters (um)
% Usage: surfaceArea = compute_surface_area (lengths, diameters, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_surface_area([3, 2], [1, 0.5])
%       compute_surface_area([3; 2], [1; 0.5])
%       compute_surface_area([3, 2], [1, 0.5], 'EachCompartment', true)
%
% Outputs:
%       surfaceArea - surface area excluding ends (cm^2)
%                   specified as a nonnegative vector
%
% Arguments:
%       lengths     - lengths (um) for each compartment (each column)
%                   Note: Each row is computed independently
%                   must be a nonnegative array
%       diameters   - diameters (um) for each compartment (each column)
%                   Note: Each row is computed independently
%                   must be a nonnegative array
%       varargin    - 'EachCompartment': whether to compute 
%                                   for each compartment separately
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Used by:
%       cd/compute_total_current.m
%       cd/m3ha_neuron_create_initial_params.m

% File History:
% 2018-11-13 Adapted from ~/m3ha/optimizer4gabab/singleneuronfitting41.m
% 

%% Hard-coded constants
CM_PER_UM = 1e-4;

%% Default values for optional arguments
eachCompartmentDefault = false;     % compute the surface area for 
                                    %   the entire cell by default

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
addRequired(iP, 'lengths', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'diameters', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'EachCompartment', eachCompartmentDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, lengths, diameters, varargin{:});
eachCompartment = iP.Results.EachCompartment;

%% Do the job
% Compute the area of the model cell in cm^2
if eachCompartment
    surfaceArea = pi .* diameters .* lengths .* (CM_PER_UM)^2;
else
    surfaceArea = sum(pi .* diameters .* lengths .* (CM_PER_UM)^2, 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
