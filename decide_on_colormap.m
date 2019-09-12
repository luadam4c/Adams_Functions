function colorMap = decide_on_colormap (colorMap, varargin)
%% Decides on the color map to use
% Usage: colorMap = decide_on_colormap (colorMap, (opt) nColors, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       decide_on_colormap([])
%       decide_on_colormap('Gray')
%       decide_on_colormap({'Red', 'Blue', 'Green'})
%       decide_on_colormap([], 4)
%       decide_on_colormap([], 4, 'ColorMapFunc', @hsv)
%
% Outputs:
%       colorMap    - color map created
%                   specified as a nColors by 3 numeric array
%
% Arguments:
%       colorMap    - color map passed in
%                   must be empty or a string/character vector
%                       or an n-by-3 numeric array
%       nColors     - (opt) number of colors
%                   must be a positive integer vector
%                   default == 64
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the create_colormap() function
%
% Requires:
%       cd/char2rgb.m
%       cd/create_colormap.m
%       cd/create_error_for_nargin.m
%       cd/match_row_count.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_traces.m
%       cd/plot_vertical_shade.m

% File History:
% 2019-08-22 Created by Adam Lu
% 

%% Hard-coded parameters
nColorsDefault = 64;

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'colorMap');

% Add optional inputs to the Input Parser
addOptional(iP, 'nColors', nColorsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, colorMap, varargin{:});
nColors = iP.Results.nColors;
% param1 = iP.Results.param1;

% Keep unmatched arguments for the create_colormap() function
otherArguments = iP.Unmatched;

%% Do the job
if isempty(colorMap)
    % Set default color map
    colorMap = create_colormap(nColors, otherArguments);
elseif ischar(colorMap)
    % Convert to a numeric array
    colorMap = char2rgb(colorMap);
elseif isstring(colorMap) || iscellstr(colorMap)
    % Convert to a numeric vectors
    % TODO: cellorarrayfun.m
    if iscell(colorMap)
        cellorarrayfun = @cellfun;
    else
        cellorarrayfun = @arrayfun;
    end
    colorMap = cellorarrayfun(@char2rgb, colorMap, 'UniformOutput', false);
    
    % Vertically concatenate them
    colorMap = vertcat(colorMap{:});
end

% Match the number of rows in the color map to nColors
colorMap = match_row_count(colorMap, nColors, 'ExpansionMethod', 'repeat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%