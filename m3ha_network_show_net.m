function [RERE, TCRE, RETC] = m3ha_network_show_net (inFolder, varargin)
%% Shows network topology for each network (each .syn file in the inFolder)
% Usage: [RERE, TCRE, RETC] = m3ha_network_show_net (inFolder, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Arguments:
%       inFolder    - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%                   must be a character array
%       varargin    - 'OutFolder' - the name of the directory where plots will be placed
%                   must be a character array
%                   default == same as inFolder
%                   - 'FirstOnly' - whether to plot the first network only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       inFolder/*.syn
%       cd/apply_iteratively.m
%       cd/set_figure_properties.m
%
% Used by:
%       cd/m3ha_launch1.m

% File History:
% 2017-10-23 Modified from /RTCl/show_net.hoc
% 2017-10-31 Now uses dir instead of dirr


%% Default values for optional arguments
firstOnlyDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help m3ha_network_show_net'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'inFolder', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutFolder', '', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'FirstOnly', firstOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, inFolder, varargin{:});
outFolder = iP.Results.OutFolder;
firstOnly = iP.Results.FirstOnly;

% Set dependent argument defaults
if isempty(outFolder)
    outFolder = inFolder;
end

%% Find all RERE*.syn files
REREfiles = dir(fullfile(inFolder, 'RERE*.syn'));
TCREfiles = dir(fullfile(inFolder, 'TCRE*.syn'));
RETCfiles = dir(fullfile(inFolder, 'RETC*.syn'));
nNetworks = length(REREfiles);          % number of networks to plot
if length(TCREfiles) ~= nNetworks || ...
    length(RETCfiles) ~= nNetworks
    error('Number of RERE files, TCRE files and RETC files must be equal!');
end

%% Change the number of networks to plot if FirstOnly flag is true
if firstOnly
    nNetworks = 1;
end

%% Create raster plot for each file
RERE = cell(nNetworks, 1);        % raw RE-RE connections
TCRE = cell(nNetworks, 1);        % raw TC-RE connections
RETC = cell(nNetworks, 1);        % raw RE-TC connections
for i = 1:nNetworks
% parfor i = 1:nNetworks
    % Extract names of files needed
    REREfilename = REREfiles(i).name;
    TCREfilename = replace(REREfilename, 'RERE', 'TCRE');
    RETCfilename = replace(REREfilename, 'RERE', 'RETC');

    % Load data
    RERE{i} = load(fullfile(inFolder, REREfilename));
    TCRE{i} = load(fullfile(inFolder, TCREfilename));
    RETC{i} = load(fullfile(inFolder, RETCfilename));

    % Set figure name
    figName = replace(REREfilename, '.syn', '.png');
    figName = replace(figName, 'RERE', 'Topology');
    figName = fullfile(outFolder, figName);

    % Plot data
    [~, fileBase, ~] = fileparts(replace(REREfilename, 'RERE_', ''));
    network_topology(RERE{i}, TCRE{i}, RETC{i}, figName, 'full', fileBase);
    network_topology(RERE{i}, TCRE{i}, RETC{i}, figName, 'part', fileBase);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function network_topology(RERE, TCRE, RETC, figName, plotMode, fileBase)
%% Plot network topology from connection data

% Prepare the figure
h = set_figure_properties('AlwaysNew', true);
hold on;

% Define color map
cm = colormap(lines);

% Count the number of colors
nColors = size(cm, 1);

% Plot TC-RE connections
for i = 1:size(TCRE, 1)
    colorThis = cm(mod(TCRE(i, 1), nColors) + 1, :);
    plot([2; 1], round([TCRE(i, 1); TCRE(i, 2)]), '-o', 'Color', colorThis);
end

% Plot RE-TC connections
for i = 1:size(RETC, 1)
    colorThis = cm(mod(RETC(i, 1), nColors) + 1, :);
    plot([1; 2], round([RETC(i, 1); RETC(i, 2)]), '-.o', 'Color', colorThis);
end

% Plot RE-RE connections
for i = 1:size(RERE, 1)
    colorThis = cm(mod(RERE(i, 1), nColors) + 1, :);
    plot([1; 0.75], round([RERE(i, 1); RERE(i, 2)]), ':', 'Color', colorThis);
end

% Restrict axes range
ymax = apply_iteratively(@max, {RERE; TCRE; RETC}) + 1;
ymin = apply_iteratively(@min, {RERE; TCRE; RETC}) - 1;
if strcmp(plotMode, 'full')
    axis([0.5, 2.5, ymin, ymax]);
elseif strcmp(plotMode, 'part')
    axis([0.5, 2.5, ymin, ymin + 11]);
end

% Set labels
ax = gca;
set(ax, 'XTick', [1, 2]);
set(ax, 'XTickLabel', {'RE', 'TC'});
ylabel('Neuron number');
title(['Network topology for ', replace(fileBase, '_', '\_')]);

% Save figure
if strcmp(plotMode, 'part')
    figName = replace(figName, '.png', '_zoomed.png');
end
saveas(h, figName, 'png');

% Close figure
% close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

ymax = max([RERE(:, 1); TCRE(:, 1); RETC(:, 1)]) + 1;
ymin = min([RERE(:, 1); TCRE(:, 1); RETC(:, 1)]) - 1;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%