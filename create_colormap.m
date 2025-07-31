function colorMap = create_colormap (varargin)
%% Returns colorMap based on the number of colors requested
% Usage: colorMap = create_colormap (nColors (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       create_colormap(1)
%       create_colormap(2)
%       create_colormap(3)
%       create_colormap(4)
%       create_colormap(7)
%       create_colormap(7, 'ColorMapFunc', @parula)
%       create_colormap(7, 'ReverseOrder', true)
%       create_colormap('ColorMapFunc', @gray, 'ReverseOrder', true)
%       create_colormap('ColorMapFunc', @gray, 'ReverseOrder', true, 'HighContrast', true)
%
% Outputs:
%       colorMap    - color map created
%                   specified as a nColors by 3 numeric array
% Arguments:
%       nColors     - (opt) number of colors
%                   must be a positive integer vector
%                   default == 64
%       varargin    - 'ColorMapFunc': color map function to use
%                       must be a function handle
%                       default == @jet
%                   - 'ReverseOrder': whether to reverse the order of the map
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'HighContrast': whether to use a high contrast map
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       ~/Downloaded_Functions/rgb.m
%       cd/locate_functionsdir.m
%
% Used by:
%       cd/decide_on_colormap.m

% File History:
% 2018-10-29 Adapted from code in run_neuron_once_4compgabab.m
% 2019-01-08 Now returns a cell array if the input has more than one elements
% 2019-05-14 If nColors == 1, make it blue
% 2019-08-04 Added colorMapFunc as an optional argument
% 2019-08-04 Made nColors an optional argument
% 2019-08-04 Added 'ReverseOrder' as an optional argument
% 2019-08-04 Added 'HighContrast' as an optional argument
% 2019-10-07 Changed default colors
% 2020-02-08 Fixed bug
% 

%% Hard-coded parameters

%% Default values for optional arguments
nColorsDefault = 64;
colorMapFuncDefault = @jet;
reverseOrderDefault = false;        % don't reverse order by default
highContrastDefault = false;        % don't use high contrast by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if exist('rgb.m', 'file') ~= 2 && ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for rgb.m, etc.
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions')); 
end

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add optional inputs to the Input Parser
addOptional(iP, 'nColors', nColorsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ColorMapFunc', colorMapFuncDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));
addParameter(iP, 'ReverseOrder', reverseOrderDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'HighContrast', highContrastDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
nColors = iP.Results.nColors;
colorMapFunc = iP.Results.ColorMapFunc;
reverseOrder = iP.Results.ReverseOrder;
highContrast = iP.Results.HighContrast;

%% Do the job
if numel(nColors) > 1
    colorMap = arrayfun(@(a) create_colormap_helper(a, colorMapFunc, ...
                            reverseOrder, highContrast), ...
                        nColors, 'UniformOutput', false);
else
    colorMap = create_colormap_helper(nColors, colorMapFunc, ...
                                     reverseOrder, highContrast);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colorMap = create_colormap_helper (nColors, colorMapFunc, ...
                                            reverseOrder, highContrast)
%% Create a single color map 
%   Note: this is an N by 3 array, where N is the number of colors

if nColors == 1
    colorMap = rgb('Blue');
elseif nColors == 2
    colorMap = [rgb('Blue'); rgb('Red')];
elseif nColors == 3
    % Color groups correspond to 3 vHold conditions
    colorMap = [rgb('Blue'); rgb('DarkGreen'); rgb('Red')];
elseif nColors == 4
    % Color groups correspond to 4 pharm conditions
    colorMap = [rgb('Black'); rgb('Blue'); ...
                rgb('Red'); rgb('Purple')];
elseif nColors == 6
    % Color groups correspond to 4 pharm conditions
    colorMap = [rgb('Blue'); rgb('Cyan'); rgb('Purple'); ...
                rgb('DarkGreen'); rgb('Red'); rgb('Orange');];
else
    % Color groups corresponding to pharm-g incr pairs
    colorMap = colorMapFunc(nColors);
end

% Reverse order if requested
if reverseOrder
    colorMap = flipud(colorMap);
end

% Apply high contrast
if highContrast
    if isequal(colorMapFunc, @gray)
        % Note: this is Mark Beenhakker's algorithm
        % Set cutoff #1
        cutOff1 = 10;
        
        % Set cutoff #2
        cutOff2 = 21;
        
        % Set cutoff #3
        cutOff3 = 29;

        % 
        colorMap(cutOff1:cutOff2-1, :) = 0.4;

        % 
        colorMap(cutOff2:cutOff3-1, :) = 0.3;

        % Make everything above cutoff #3 black
        colorMap(cutOff3:end, :) = 0;
    else
        % TODO
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

colorMap = colormap(colorMapFunc(nColors));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
