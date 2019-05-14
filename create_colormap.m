function colorMap = create_colormap (nColors, varargin);
%% Returns colorMap based on the number of colors requested
% Usage: colorMap = create_colormap (nColors, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       colorMap     - TODO: Description of colorMap
%                   specified as a TODO
% Arguments:
%       nColors     - number of colors
%                   must be a positive integer vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       ~/Downloaded_Functions/rgb.m
%       cd/locate_functionsdir.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/plot_traces.m   

% File History:
% 2018-10-29 Adapted from code in run_neuron_once_4compgabab.m
% 2019-01-08 Now returns a cell array if the input has more than one elements
% 2019-05-14 If nColors == 1, make it blue
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for rgb.m, etc.
    addpath(fullfile(functionsDirectory, 'Downloaded_Functions')); 
end

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
addRequired(iP, 'nColors', ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, nColors, varargin{:});

%% Do the job
if numel(nColors) > 1
    colorMap = arrayfun(@create_colormap_helper, nColors, ...
                        'UniformOutput', false);
else
    colorMap = create_colormap_helper(nColors);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colorMap = create_colormap_helper (nColors)
%% Create a single color map 
%   Note: this is an N by 3 array, where N is the number of colors

if nColors == 1
    colorMap = rgb('Blue');    
elseif nColors == 3
    % Color groups correspond to 3 vHold conditions
    colorMap = [rgb('Blue'); rgb('Red'); rgb('Purple')];
elseif nColors == 4
    % Color groups correspond to 4 pharm conditions
    colorMap = [rgb('Black'); rgb('Blue'); ...
                rgb('Red'); rgb('Purple')];
else
    % Color groups corresponding to pharm-g incr pairs
    colorMap = colormap(jet(nColors));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
