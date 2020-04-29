function [radiusSoma, radiusDend, lengthDend] = decide_on_geom_params (varargin)
%% Standardizes geometric parameters
% Usage: [radiusSoma, radiusDend, lengthDend] = decide_on_geom_params (geomParams (opt), radiusSoma (opt), radiusDend (opt), lengthDend (opt), varargin)
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
%       geomParams  - (opt) geometric parameters
%                   must be a scalar struct with fields:
%                       radiusSoma or diamSoma
%                       radiusDendrite or radiusDend or diamDendrite or diamDend
%                       lengthDendrite or lengthDend or LDend
%                   default == struct.empty
%       radiusSoma  - (opt) radius of soma
%                   must be a nonegative scalar
%       radiusDend  - (opt) radius of dendrite
%                   must be a nonegative scalar
%       lengthDend  - (opt) length of dendrite
%                   must be a nonegative scalar
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/first_matching_field.m
%
% Used by:
%       cd/plot_ball_stick.m

% File History:
% 2019-12-26 Now uses first_matching_field.m
% 2019-12-31 Moved from by plot_ball_stick.m
% 2020-01-31 Added 'soma_radius', 'soma_diam', 'dend1_diam', 'dend1_L'

%% Hard-coded parameters
radiusSomaStr = {'radiusSoma', 'soma_radius'};
diamSomaStr = {'diamSoma', 'soma_diam'};
radiusDendriteStr = {'radiusDendrite', 'radiusDend'};
diamDendriteStr = {'diamDendrite', 'diamDend', 'dend1_diam'};
lengthDendriteStr = {'lengthDendrite', 'lengthDend', 'LDend'};
lengthDendrite1Str = {'dend1_L'};

%% Default values for optional arguments
geomParamsDefault = struct.empty;
radiusSomaDefault = [];
radiusDendDefault = [];
lengthDendDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add optional inputs to the Input Parser
addOptional(iP, 'geomParams', geomParamsDefault, ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addOptional(iP, 'radiusSoma', radiusSomaDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));
addOptional(iP, 'radiusDend', radiusDendDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));
addOptional(iP, 'lengthDend', lengthDendDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));

% Read from the Input Parser
parse(iP, varargin{:});
geomParams = iP.Results.geomParams;
radiusSoma = iP.Results.radiusSoma;
radiusDend = iP.Results.radiusDend;
lengthDend = iP.Results.lengthDend;

%% Do the job
% Radius of soma
if isempty(radiusSoma)
    if ~isempty(geomParams)
        % Look for a somatic radius
        radiusSoma = first_matching_field(geomParams, radiusSomaStr);

        % If not found, look for a dendritic diameter
        if isempty(radiusSoma)
            % Look for a somatic diameter
            diamSoma = first_matching_field(geomParams, diamSomaStr);

            % Compute the somatic radius
            if ~isempty(diamSoma)
                radiusSoma = diamSoma ./ 2;
            else
                error('No soma radius provided in geomParams!');
            end
        end
    else
        error('No soma radius passed in!');
    end
end

% Radius of dendrite
if isempty(radiusDend)
    if ~isempty(geomParams)
        % Look for a dendritic radius
        radiusDend = first_matching_field(geomParams, radiusDendriteStr);

        % If not found, look for a dendritic diameter
        if isempty(radiusDend)
            % Look for a dendritic diameter
            diamDend = first_matching_field(geomParams, diamDendriteStr);

            % Compute the dendritic radius
            if ~isempty(diamDend)
                radiusDend = diamDend ./ 2;
            else
                error('No dendrite radius provided in geomParams!');
            end
        end
    else
        error('No dendrite radius passed in!');
    end
end

% Length of dendrite
if isempty(lengthDend)
    if ~isempty(geomParams)
        % Look for a dendritic length
        lengthDend = first_matching_field(geomParams, lengthDendriteStr);

        % If not found, return error
        if isempty(lengthDend)
            % Look for a dendrite 1 length
            %   Note: this assumes there are only two dendrites of equal length
            lengthDend1 = first_matching_field(geomParams, lengthDendrite1Str);

            if ~isempty(lengthDend1)
                lengthDend = lengthDend1 .* 2;
            else
                error('No dendrite length provided in geomParams!');
            end
        end
    else
        error('No dendrite length passed in!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%