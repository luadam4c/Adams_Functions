function h = plot_tuning_map (pValues, readout, varargin)
%% Plot a 2-dimensional tuning map
% Usage: h = plot_tuning_map (pValues, readout, varargin)
% Outputs:
%       h           - figure handle for the created figure
%                   specified as a figure handle
% Arguments:
%       pValues     - 2 vectors of parameter values
%                   must be a 2-element cell array of numeric vectors
%       readout     - a readout matrix corresponding to the parameter-pairs
%                   must be a 2-dimensional numeric array with dimensions 
%                   given by the lengths of elements of pValues
%       varargin    - 'PisLog': whether parameter values are to be plotted 
%                               log-scaled
%                   must be a 2-element binary array
%                   default == [false, false];
%                   - 'PLabels': labels for parameters
%                   must be a 2-element cell array of string scalars or 
%                       character vectors
%                   default == {'Parameter1', 'Parameter2'}
%                   - 'ReadoutLabel': label for readout matrix
%                   must be a string scalar or a character vector
%                   default == 'Readout'
%                   - 'XLimits': limits of x axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'YLimits': limits of y axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'CLimits': limits of color axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                                   e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in 
%                       saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%
% Used by:    
%       cd/m3ha_network_tuning_maps.m
%       /media/adamX/RTCl/tuning_maps.m

% File History:
% 2017-05-05 Created by Adam Lu
% 2017-05-09 Added 'FigTypes' as a parameter-value pair argument
% 2017-12-18 Changed tabs to spaces
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 

%% Default values for optional arguments
pislogDefault = [false, false];
pLabelsDefault = {'Parameter1', 'Parameter2'};
readoutLabelDefault = 'Readout';
xlimitsDefault = [];
ylimitsDefault = [];
climitsDefault = [];
figNameDefault = '';
figTypesDefault = 'png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an Input Parser
addRequired(iP, 'pValues', ...              % 2 vectors of parameter values
    @(x) assert(iscell(x) && numel(x) == 2 ...
        && isnumeric(x{1}) && isvector(x{1}) ...
        && isnumeric(x{2}) && isvector(x{1}), ...
        'pValues must be a 2-element cell array of numeric vectors!'));
addRequired(iP, 'readout', ...              % a readout matrix
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PisLog', pislogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, ...
                                {'binary', 'vector', 'numel', 2}));
addParameter(iP, 'PLabels', pLabelsDefault, ...
    @(x) assert(iscell(x) && numel(x) == 2 ...
        && (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        ['PLabels must be a 2-element cell array of ', ...
            'string scalars or character vectors!']));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLimits', xlimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                            {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'YLimits', ylimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                            {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'CLimits', climitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                            {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, pValues, readout, varargin{:});
pIsLog = iP.Results.PisLog;
pLabels = iP.Results.PLabels;
readoutLabel = iP.Results.ReadoutLabel;
xlimits = iP.Results.XLimits;
ylimits = iP.Results.YLimits;
climits = iP.Results.CLimits;
figName = iP.Results.FigName;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Check relationships between arguments
if size(readout) ~= cellfun(@length, pValues)
    error(['Readout matrix must have dimensions given ', ...
            'by the lengths of elements of pValues!']);
end

%% Prepare for tuning curve
% Decide on the figure to plot on
if ~isempty(figName)
    % Create and clear figure if figName provided
    h = figure(floor(rand()*10^4)+1);
    clf(h);
else
    % Get the current figure
    h = gcf;
end

%% Plot tuning map
% Plot readout values against parameter values
surf(pValues{1}, pValues{2}, readout');
if pIsLog(1)
    set(gca, 'XScale', 'log');
elseif pIsLog(2)
    set(gca, 'YScale', 'log');
end
set(gca, 'CLimits', climits);
colorbar;

% Restrict x axis if xlimits provided
if ~isempty(xlimits)
    xlim(xlimits);
end

% Restrict y axis if ylimits provided
if ~isempty(ylimits)
    ylim(ylimits);
end

% Set title and axes labels; allow to be suppressed
if ~isequal(pLabels{1}, 'suppress')
    xlabel(pLabels{1});
end
if ~isequal(pLabels{2}, 'suppress')
    ylabel(pLabels{2});
end
if ~isequal(pLabels{1}, 'suppress') && ...
    ~isequal(pLabels{2}, 'suppress') && ~isequal(readoutLabel, 'suppress')
    title(strrep(readoutLabel, '_', '\_'));
%    title(strrep([readoutLabel, ' vs. ', pLabels{1}, ...
%                   ' and ', pLabels{2}], '_', '\_'));
end

%% Post-plotting
% Save figure if figName provided
if ~isempty(figName)
    save_all_figtypes(h, figName, figtypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%