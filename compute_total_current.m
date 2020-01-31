function [totalCurrent, compCurrents] = ...
                compute_total_current (compCurrentDensities, varargin)
%% Computes the total current across all compartments from current densities, lengths and diameters of each compartment
% Usage: [totalCurrent, compCurrents] = ...
%               compute_total_current (compCurrentDensities, lengths, diameters, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_total_current([10, 30], [3, 2], [1, 0.5])
%
% Outputs:
%       totalCurrent    - total current (nA) across all compartments
%                       specified as a numeric scalar
%
% Arguments:
%       compCurrentDensities    - current densities (mA/cm2) 
%                                       in each compartment (each column)
%                               Note: Each row is computed independently
%                               must be a numeric vector
%       lengths         - lengths (um) of each compartment (each column)
%                       Note: Each row is computed independently
%                       must be a numeric vector
%       diameters       - diameters (um) of each compartment (each column)
%                       Note: Each row is computed independently
%                       must be a numeric vector
%       varargin    - 'GeomParams': geometric parameters
%                   must be a scalar struct with fields:
%                       radiusSoma
%                       radiusDendrite
%                       lengthDendrite
%                   default == struct.empty
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/compute_surface_area.m
%
% Used by:
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_plot_simulated_traces.m

% File History:
% 2019-12-31 Created by Adam Lu
% 

%% Hard-coded parameters
NA_PER_MA = 1e6;

%% Default values for optional arguments
lengthsDefault = [];
diametersDefault = [];
geomParamsDefault = struct.empty;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'compCurrentDensities', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'lengths', lengthsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addOptional(iP, 'diameters', diametersDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'GeomParams', geomParamsDefault, ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Read from the Input Parser
parse(iP, compCurrentDensities, varargin{:});
lengths = iP.Results.lengths;
diameters = iP.Results.diameters;
geomParams = iP.Results.GeomParams;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Preparation
% TODO: m3ha_decide_on_geom_params.m
% Read from geomParams if needed
if isempty(lengths) || isempty(diameters)
    [radiusSoma, radiusDend, lengthDend] = decide_on_geom_params(geomParams);
end

% Decide on the lengths
if isempty(lengths)
    lengths = [radiusSoma * 2, lengthDend / 2, lengthDend / 2];
end

% Decide on the diameters
if isempty(diameters)
    diameters = [radiusSoma, radiusDend, radiusDend] * 2;
end

%% Do the job
% Compute the surface area (cm^2) of each compartment separately
surfaceArea = compute_surface_area(lengths, diameters, 'EachCompartment', true);

% Compute the component currents in nA
compCurrents = compCurrentDensities .* surfaceArea .* NA_PER_MA;

% Compute the total current in nA
totalCurrent = sum(compCurrents, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%