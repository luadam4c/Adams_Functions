function m3ha_histograms_across_cells (identifier, varargin)
%% Plots histograms across cells for single neuron fitting
% Usage: m3ha_histograms_across_cells (identifier, varargin)
% Explanation:
%       Plots error histograms
% Arguments:
%   identifier  - must be character array/string
%   varargin    - 'ColorCode': Color code cells w/ legend
%               must be numeric/logical 1 (true) or 0 (false)
%                   0    - plots without color
%                   1    - plots with color
%               default == 0 XXX: Should default be color or not?
%               'IncludeInit': Plot initial error
%               must be numeric/logical 1 (true) or 0 (false)
%                   0   - plots only optimized error
%                   1   - plots both initial and optimized error
%               default == 0
%               if IncludeInit == 1, ColorCode is unused
%
% File History:
% 2017-08-23 - BT - Created by Brian
% 2017-08-25 - AL - Changed name error_histogram.m -> histograms_across_cells.m
% 2017-08-25 - BT - snfVersion -> identifier, modified filename structuring, extracted fields
% 2017-08-28 - BT - histg for color option, converted csv read and histogram flow for parfor
% 2017-08-30 - BT - Changed labels, optional -> parameters
% 2017-09-01 - BT - Removed length(unique(x)) == 1 fields from plotting
% 2018-01-24 - AL - Removed addpath - Is it needed?
% XXX: Error from bar.m in histg.m: XData cannot contain duplicate values. Arrays for which length(unique(x)) == 1 do not all have this problem. Only at certain values (corrD @ 7.9540, and having checked 7.9530-7.9540, 7.9538 and 7.9530) do this occur. Tested manually using histg.

%% Parameters used in the body of the function
csvidentifier = 'errors_and_params_log_manual_concise.csv';

%% Default values for optional arguments
ColorCodeDefault = 0;        % default ColorCode: do not color code cells w/ legend
IncludeInitDefault  = 0;    % default IncludeInit: do not include initial error

%% Directory (ies) for placing outputs
outfolder_suf = 'histograms_across_cells';
outfolder = ['/media/adamX/m3ha/optimizer4gabab/' identifier '_' outfolder_suf '/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with directories
% Find home directory across servers
if exist('/media/adamX/m3ha/', 'dir') == 7
    homeDirectory = '/media/adamX/m3ha/';
elseif exist('/scratch/al4ng/m3ha/', 'dir') == 7
    homeDirectory = '/scratch/al4ng/m3ha/';
else
    error('Valid homeDirectory does not exist!');
end

%% Deal with arguments
% Check number of required arguments
if nargin < 1
   error(['Not enough input arguments, ', ...
            'type ''help histogram_across_cells'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'histogram_across_cells';
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'identifier', ...               % identifier must be char array/string
    @(x) validateattributes(x, {'char'}, {'nonempty' 'row'}));

% Add optional inputs to the Input Parser
addParameter(iP, 'ColorCode', ColorCodeDefault, ... % ColorCode is 0 or 1
    @(x) assert((x == 0 || x == 1) && all(size(x) == 1)));

addParameter(iP, 'IncludeInit', IncludeInitDefault, ...  % IncludeInit is 0 or 1
    @(x) assert((x == 0 || x == 1) && all(size(x) == 1)));

% Read from the Input Parser
parse(iP, identifier, varargin{:});
identifier = iP.Results.identifier;
ColorCode = iP.Results.ColorCode;
IncludeInit = iP.Results.IncludeInit;

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

% Check if needed output directories exist
if exist(fullfile(outfolder), 'dir') ~= 7
    newdirectory = fullfile(outfolder);
    mkdir(newdirectory);
    fprintf('New directory made: %s\n\n', newdirectory);
end

% Find folders and csv files matching identifier
dirStruct = dir([identifier '*']);  % matching identifier
dirCell = struct2cell(dirStruct);   % convert to cell array
dirAll = dirCell(1,:);  % extract names in dir
dirFolders = dirAll(cellfun(@(x) isdir(x), dirAll));    % get folders in dir

potIndex = strfind(dirFolders, outfolder_suf);    % remove possible duplicate created by outfolder
dupIndex = find(not(cellfun('isempty', potIndex)));
dirFolders(dupIndex) = [];

indexID = strfind(dirFolders{1}, '_');  % position of underscore to get ID
dirFoldersID = cellfun(@(x) x(indexID(1)+1:indexID(2)-1), dirFolders, 'UniformOutput', false);    % get unique ID
dirExcel = cellfun(@(x) [identifier '_' x '_' csvidentifier], dirFoldersID, 'UniformOutput', false); % filenames of csv files

fields = readtable(fullfile(dirFolders{1}, dirExcel{1}));    % get fields of first file
fields = fields.Properties.VariableNames;    % extract variable names

% Read files and plot
allData = [];
origData = [];
if IncludeInit == 0 % initial error not included
    parfor l = 1:size(dirFoldersID,2)
        allData(l,:) = csvread(fullfile(dirFolders{l}, dirExcel{l}), 2, 0); % reads row of optimized
        tempfields = readtable(fullfile(dirFolders{l}, dirExcel{l}));    % compare fields of current file to first
        tempfields = tempfields.Properties.VariableNames;
        if ~isempty(setdiff(fields, tempfields)) || ~isempty(setdiff(tempfields, fields))
            disp([dirFoldersID ' has missing/extraneous fields!']);
        end
    end
    
    parfor h = 1:size(allData,2)
        f = figure(h);
        f.Visible = 'off';
        if ColorCode == 0    % plot without color
            histogram(allData(:,h));    % plot error for one field
        elseif ColorCode == 1    % plot with color and legend
            classes = 1:size(dirFolders,2);    % unique class for each
            histg(allData(:,h)', classes);    % histg for one field, transposed data for histg
            legend(dirFoldersID, 'location', 'eastoutside');    % label only for color
        end
        title(['Optimized ' fields{h}]);   % labels and saving
        xlabel(fields{h});
        ylabel(['Frequency of ' fields{h}]);
        figname = fullfile(outfolder, [identifier '_' fields{h} '.png']);
        saveas(f, figname);
    end
elseif IncludeInit == 1 % initial error included
    parfor l = 1:size(dirFoldersID,2)
        origData(l,:) = csvread(fullfile(dirFolders{l}, dirExcel{l}), 1, 0, [1 0 1 size(fields,2)-1]); % reads row of original
        optData(l,:) = csvread(fullfile(dirFolders{l}, dirExcel{l}), 2, 0); % reads row of optimized
        tempfields = readtable(fullfile(dirFolders{l}, dirExcel{l}));    % compare fields of current file to first
        tempfields = tempfields.Properties.VariableNames;
        if ~isempty(setdiff(fields, tempfields)) || ~isempty(setdiff(tempfields, fields))
            disp([dirFoldersID ' has missing/extraneous fields!']);
        end
    end
    allData = [origData; optData];    % combine data
    opt.group_names(1:size(origData,1),1) = {'Original'};    % group names cell array for Original and Optimized groups
    opt.group_names(size(origData,1)+1:size(origData,1)+size(optData,1),1) = {'Optimized'};
    classes = [ones(1,size(origData,1)) size(allData,1)*ones(1,size(optData,1))];    % class for original then optimized data, size(allData,1) determines the color of optimized
    parfor h = 1:size(allData,2)
        if size(unique(allData(:,h)),1) ~= 1    % prevents errors from corrD error 
            f = figure(h);
            f.Visible = 'off';
            z = histg(allData(:,h)', classes, opt);    % histg for one field, transposed data for histg
            legend([z(1) z(end)], {'Original', 'Optimized'}, 'location', 'eastoutside');    % label only for color, [z(1) z(end)] selects bars to color
            title(['Original and optimized ' fields{h}]);   % labels and saving
            xlabel(fields{h});
            ylabel(['Frequency of ' fields{h}]);
            figname = fullfile(outfolder, [identifier '_' fields{h} '_both.png']);
            saveas(f, figname);
        end
    end
    for g = 1:size(origData,1)  % groups together data by cell
        classes(g) = g;
        classes(g+size(origData,1)) = g;
        opt.group_names{g} = dirFoldersID{g};
        opt.group_names{g+size(origData,1)} = dirFoldersID{g};
    end
    parfor h = 1:size(allData,2)
        if size(unique(allData(:,h)),1) ~= 1
            f = figure(h);
            histg(allData(:,h)', classes, opt);
            legend(dirFoldersID, 'location', 'eastoutside');
            title(['Original and optimized ' fields{h}]);   % labels and saving
            xlabel(fields{h});
            ylabel(['Frequency of ' fields{h}]);
            figname = fullfile(outfolder, [identifier '_' fields{h} '_both_bycell.png']);
            saveas(f, figname);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:
TODO: Place older versions of the code that you want to save here, 
        in case you need it back in the future
    parfor l_iter = 1:size(dirFoldersID,2)
        allData(l:l+1,:) = csvread(fullfile(dirFolders{l}, dirExcel{l}), 1, 0); % reads row of optimized
        tempfields = readtable(fullfile(dirFolders{l}, dirExcel{l}));    % compare fields of current file to first
        tempfields = tempfields.Properties.VariableNames;
        if ~isempty(setdiff(fields, tempfields)) || ~isempty(setdiff(tempfields, fields))
            disp([dirFoldersID ' has missing/extraneous fields!']);
        end
    end
        f = figure(h);
        f.Visible = 'off';
        classes = (-1).^(0:size(allData,1));    % alternating class for original and optimized
        histg(allData(:,h)', classes);    % histg for one field, transposed data for histg
        legend({'Original', 'Optimized'}, 'location', 'eastoutside');    % label only for color
        title(['Original and optimized ' fields{h}]);   % labels and saving
        xlabel(fields{h});
        ylabel(['Frequency of ' fields{h}]);
        figname = fullfile(outfolder, [identifier '_' fields{h} '_both.png']);
        saveas(h, figname);
%}
