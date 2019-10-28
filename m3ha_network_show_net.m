function [RERE, TCRE, RETC] = m3ha_network_show_net (infolder, varargin)
%% Shows network topology for each network (each .syn file in the infolder)
% USAGE: [RERE, TCRE, RETC] = m3ha_network_show_net (infolder, varargin)
%
% Arguments:
%       infolder    - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%                   must be a character array
%       varargin    - 'OutFolder' - the name of the directory where plots will be placed
%                   must be a character array
%                   default == same as infolder
%                   - 'FirstOnly' - whether to plot the first network only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%        infolder/*.syn
% Used by:
%        cd/m3ha_launch1.m
%
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
addRequired(iP, 'infolder', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutFolder', '', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'FirstOnly', firstOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, infolder, varargin{:});
outfolder = iP.Results.OutFolder;
firstonly = iP.Results.FirstOnly;

% Set dependent argument defaults
if isempty(outfolder)
    outfolder = infolder;
end

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end

%% Find all RERE*.syn files
REREfiles = dir(fullfile(infolder, 'RERE*.syn'));
TCREfiles = dir(fullfile(infolder, 'TCRE*.syn'));
RETCfiles = dir(fullfile(infolder, 'RETC*.syn'));
nNetworks = length(REREfiles);          % number of networks to plot
if length(TCREfiles) ~= nNetworks || ...
    length(RETCfiles) ~= nNetworks
    error('Number of RERE files, TCRE files and RETC files must be equal!');
end

%% Change the number of networks to plot if FirstOnly flag is true
if firstonly
    nNetworks = 1;
end

%% Create raster plot for each file
RERE = cell(nNetworks, 1);        % raw RE-RE connections
TCRE = cell(nNetworks, 1);        % raw TC-RE connections
RETC = cell(nNetworks, 1);        % raw RE-TC connections
parfor i = 1:nNetworks
    % Extract names of files needed
    REREfilename = REREfiles(i).name;
    TCREfilename = strrep(REREfilename, 'RERE', 'TCRE');
    RETCfilename = strrep(REREfilename, 'RERE', 'RETC');

    % Load data
    RERE{i} = load(fullfile(infolder, REREfilename));
    TCRE{i} = load(fullfile(infolder, TCREfilename));
    RETC{i} = load(fullfile(infolder, RETCfilename));

    % Set figure name
    figname = strrep(REREfilename, '.syn', '.png');
    figname = strrep(figname, 'RERE', 'Topology');
    figname = fullfile(outfolder, figname);

    % Plot data
    [~, filebase, ~] = fileparts(strrep(REREfilename, 'RERE_', ''));
    network_topology(RERE{i}, TCRE{i}, RETC{i}, figname, 'full', filebase);
    network_topology(RERE{i}, TCRE{i}, RETC{i}, figname, 'part', filebase);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function network_topology(RERE, TCRE, RETC, figname, plotmode, filebase)
%% Plot network topology from connection data

% Prepare the figure
h = figure(floor(rand()*10^4));
clf(h);
hold on;

% Define color map
cm = colormap(lines);
ncolors = size(cm, 1);

% Plot TC-RE connections
for i = 1:length(TCRE)
    color = cm(mod(TCRE(i, 1), ncolors) + 1, :);
    if strcmp(plotmode, 'full')
        plot([2; 1], round([TCRE(i, 1); TCRE(i, 2)]), 'r-o', 'Color', color);
    elseif strcmp(plotmode, 'part')
        plot([2; 1], round([TCRE(i, 1); TCRE(i, 2)]), 'r-o', 'Color', color);
    end
end

% Plot RE-TC connections
for i = 1:length(RETC)
    color = cm(mod(RETC(i, 1), ncolors) + 1, :);
    if strcmp(plotmode, 'full')
        plot([1; 2], round([RETC(i, 1); RETC(i, 2)]), 'g-.o', 'Color', color);
    elseif strcmp(plotmode, 'part')
        plot([1; 2], round([RETC(i, 1); RETC(i, 2)]), 'g-.o', 'Color', color);
    end
end

% Plot RE-RE connections
for i = 1:length(RERE)
    color = cm(mod(RERE(i, 1), ncolors) + 1, :);
    if strcmp(plotmode, 'full')
        plot([1; 0.75], round([RERE(i, 1); RERE(i, 2)]), 'b:', 'Color', color);
    elseif strcmp(plotmode, 'part')
        plot([1; 0.75], round([RERE(i, 1); RERE(i, 2)]), 'b:', 'Color', color);
    end
end

% Restrict axes range
ymax = max([RERE(:, 1); TCRE(:, 1); RETC(:, 1)]) + 1;
ymin = min([RERE(:, 1); TCRE(:, 1); RETC(:, 1)]) - 1;
if strcmp(plotmode, 'full')
    axis([0.5 2.5 ymin ymax]);
elseif strcmp(plotmode, 'part')
    axis([0.5 2.5 ymin ymin+11]);
end

% Set labels
ax = gca;
set(ax, 'XTick', [1, 2]);
set(ax, 'XTickLabel', {'RE', 'TC'});
ylabel('Neuron number');
title(['Network topology for ', strrep(filebase, '_', '\_')]);

% Save figure
if strcmp(plotmode, 'part')
    figname = strrep(figname, '.png', '_zoomed.png');
end
saveas(h, figname, 'png');

% Close figure
% close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%