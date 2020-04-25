function colorMap = decide_on_colormap (colorMap, varargin)
%% Decides on the color map to use
% Usage: colorMap = decide_on_colormap (colorMap, (opt) nColors, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       decide_on_colormap([])
%       decide_on_colormap('Gray')
%       decide_on_colormap({'Red', 'Green'}, 2)
%       decide_on_colormap([], 4)
%       decide_on_colormap(@hsv, 4)
%       decide_on_colormap({'Red', 'Blue', 'Green'})
%       decide_on_colormap({'Red', 'Blue', 'Green'}, 'ForceCellOutput', true)
%       decide_on_colormap({'Red', 'Blue', 'Green'}, 'OriginalNColors', false, 'ForceCellOutput', true)
%       decide_on_colormap({'Red', 'Blue', 'Green'}, 'FadePercentage', 50)
%
% Outputs:
%       colorMap    - color map created
%                   specified as a nColors by 3 numeric array
%
% Arguments:
%       colorMap    - color map passed in
%                   must be empty or a string/character vector
%                       or an n-by-3 numeric array
%                       or a cell array of n-by-3 numeric arrays
%                       or a function handle
%       nColors     - (opt) number of colors
%                   must be a positive integer vector
%                   default == 64
%       varargin    - 'ForceCellOutput': whether to force output as 
%                                           a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OriginalNColors': whether to keep 
%                                       original number of colors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == if forceCellOutput is true, true;
%                               otherwise, expands to nColors
%                   - 'DarkPercentage': darking out percentage
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == none
%                   - 'FadePercentage': fading out percentage
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == none
%                   - Any other parameter-value pair for 
%                       the create_colormap() function
%
% Requires:
%       cd/array_fun.m
%       cd/char2rgb.m
%       cd/create_colormap.m
%       cd/create_error_for_nargin.m
%       cd/match_row_count.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_figure07.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/plot_bar.m
%       cd/plot_chevron.m
%       cd/plot_grouped_histogram.m
%       cd/plot_grouped_jitter.m
%       cd/plot_raster.m
%       cd/plot_selected.m
%       cd/plot_spectrogram.m
%       cd/plot_spike_density_multiunit.m
%       cd/plot_traces.m
%       cd/plot_vertical_shade.m
%       cd/plot_violin.m

% File History:
% 2019-08-22 Created by Adam Lu
% 2019-10-12 Now allows colorMap to be a cell array
% 2019-12-18 Now allows colorMap to be a function handle
% 2020-04-15 Added 'ForceCellOutput' as an optional argument
% 2020-04-20 Added 'OriginalNColors' as an optional argument
% 2020-04-20 Added 'DarkPercentage' as an optional argument
% 2020-04-20 Added 'FadePercentage' as an optional argument

%% Hard-coded constants
WHITE = [1, 1, 1];

%% Hard-coded parameters
defaultNColors = 64;

%% Default values for optional arguments
nColorsDefault = [];                % set later
forceCellOutputDefault = false;     % don't force as cell array by default
originalNColorsDefault = [];        % set later
darkPercentageDefault = [];         % don't dark out by default
fadePercentageDefault = [];         % don't fade out by default

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
    @(x) validateattributes(x, {'numeric'}, {'vector', 'nonnegative', 'integer'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OriginalNColors', originalNColorsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'DarkPercentage', darkPercentageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));
addParameter(iP, 'FadePercentage', fadePercentageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));

% Read from the Input Parser
parse(iP, colorMap, varargin{:});
nColors = iP.Results.nColors;
forceCellOutput = iP.Results.ForceCellOutput;
originalNColors = iP.Results.OriginalNColors;
darkPercentage = iP.Results.DarkPercentage;
fadePercentage = iP.Results.FadePercentage;

% Keep unmatched arguments for the create_colormap() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default originalNColors
if isempty(originalNColors)
    if forceCellOutput
        originalNColors = true;
    else
        originalNColors = false;
    end
end

% Set default nColors
if isempty(nColors)
    if originalNColors && ~isempty(colorMap)
        if isstring(colorMap) || iscellstr(colorMap)
            nColors = numel(colorMap);
        elseif ischar(colorMap)
            nColors = 1;
        else
            nColors = size(colorMap, 1);
        end
    else
        nColors = defaultNColors;
    end
end

%% Do the job
if isempty(colorMap)
    % Set default color map
    colorMap = create_colormap(nColors, otherArguments);
elseif isa(colorMap, 'function_handle')
    % Set default color map
    colorMap = create_colormap(nColors, 'ColorMapFunc', colorMap, ...
                                otherArguments);
elseif ischar(colorMap)
    % Convert to a numeric array
    colorMap = char2rgb(colorMap);
elseif isstring(colorMap) || iscellstr(colorMap)
    % Convert to a numeric vectors
    colorMap = array_fun(@char2rgb, colorMap, 'UniformOutput', false);
    
    % Vertically concatenate them
    colorMap = vertcat(colorMap{:});
end

% Dark out if requested
if ~isempty(darkPercentage)
    colorMap = colorMap .* (darkPercentage / 100);
end

% Fade out if requested
if ~isempty(fadePercentage)
    colorMap = WHITE - (WHITE - colorMap) * (fadePercentage / 100);
end

% Match the number of rows in the color map to nColors
if isscalar(nColors)
    colorMap = match_row_count(colorMap, nColors, 'ExpansionMethod', 'repeat');
end

% Force as a cell array if requested
if forceCellOutput
    rowColors = transpose(1:size(colorMap, 1));
    colorMap = arrayfun(@(i) colorMap(i, :), rowColors, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%