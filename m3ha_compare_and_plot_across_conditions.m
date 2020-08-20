function m3ha_compare_and_plot_across_conditions (dirIdentifier, varargin)
%% Plot activation/inactivation and I-V curves across conditions
% Usage: m3ha_compare_and_plot_across_conditions (dirIdentifier, varargin)
%
% Arguments: TODO
%
% Requires:
%       cd/apply_over_cells.m
%       cd/m3ha_compare_neuronparams2.m

% File History:
% 2017-09-15 Adapted from m3ha_compare_and_plot_across_IC2.m
% 2017-09-20 BT - Find pFileNames in dir by suffixes and 
%                   load into m3ha_compare_neuronparams2
% 2018-08-17 AL - Improved code legibility
% 2018-08-17 AL - Made identifier an optional arguement

%% Hard-coded parameters
subDirPattern = '^\d{8}T\d{4}_\w\d{6}$';
cellNamePattern = '_\w\d{6}';

%% Default values for optional arguments
suffixesDefault = {'_bef'; '_aft'};
outFolderDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required arguments to the Input Parser
addRequired(iP, 'dirIdentifier', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Suffixes', suffixesDefault, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, dirIdentifier, varargin{:});
suffixes = iP.Results.Suffixes;
outFolder = iP.Results.OutFolder;

%% Preparation
% Decide on the home directory
% TODO: Make this an optional argument
homeDirectory = pwd;

% Decide on the input directory pattern
if isempty(dirIdentifier)
    inputDirectoryPattern = 'results';
else
    inputDirectoryPattern = dirIdentifier;
end

% Decide on the .p file regexp patterns
pfilePatterns = cellfun(@(x) [x, '\.p$'], suffixes, 'UniformOutput', false);

% Decide on the output directory
if isempty(outFolder)
    outFolder = fullfile(pwd, 'results', ...
                            ['compare_neuronparams_', dirIdentifier]);
end

%% Do the job
% List all subDirs/files within the home directory
listingHome = dir(homeDirectory);

% Extract subdirectories only
inputCands = listingHome([listingHome.isdir]);

% Extract subdirectories matching the input directory pattern only
inputDirs = ...
    inputCands(cellfun(@(x) any(regexp(x, inputDirectoryPattern)), ...
                        {inputCands.name}));

% Find all subdirectories of interest
subDirNames = {};
subDirParents = {};
for iInputDir = 1:numel(inputDirs)
    % Get the path of the current input directory
    thisInputDir = fullfile(inputDirs(iInputDir).folder, ...
                            inputDirs(iInputDir).name);

    % Get all subDirs/files within the input directories of that pattern
    listing = dir(thisInputDir);

    % Extract subdirectories only
    subDirs = listing([listing.isdir]);

    % Extract subdirectories matching the pattern only
    subDirsOfInterest = subDirs(cellfun(@(x) any(regexp(x, subDirPattern)), ...
                                        {subDirs.name}));

    % Get all subdirectory names and folders
    subDirNames = [subDirNames, {subDirsOfInterest.name}];
    subDirParents = [subDirParents, {subDirsOfInterest.folder}];
end

% Count number of subdirectories
nSubDirs = numel(subDirNames);

% Loop through all subdirectories of interest
parfor iSubDirs = 1:nSubDirs
%for iSubDirs = 1:nSubDirs
    % Get the full path to the current subdirectory
    subDirPath = fullfile(subDirParents{iSubDirs}, subDirNames{iSubDirs});

    % Get all files in the current subdirectory
    subFiles = dir(subDirPath);
    subFileNames = {subFiles.name};

    % Find all files in the current subdirectory
    %   matching the pfile patterns of interest
    pFileNamesAll = ...
        cellfun(@(y) subFileNames(cellfun(@(x) any(regexp(x, y)), ...
                                          subFileNames)), ...
                pfilePatterns, 'UniformOutput', false);

    % Put all matching pFile names together
    pFileNames = apply_over_cells(@union, pFileNamesAll, 'OptArg', 'stable');

    % Count the number of files
    nPFiles = numel(pFileNames);

    % Only do anything if pFileNames exist
    if ~isempty(pFileNames)
        % Preallocate paramNames, paramValues and suffixes
        paramNames = cell(nPFiles, 1);
        paramValues = cell(nPFiles, 1);
        suffixesActual = cell(nPFiles, 1);

        % Import p file data in turn
        for iPFile = 1:nPFiles
            % Get the current pfile name
            pFileName = pFileNames{iPFile};

            % Extract the .p file base
            [~, pFileBase, ~] = fileparts(pFileName);

            % Get the cell name string including preceding '_'
            [idxStart, idxEnd] = regexp(pFileBase, cellNamePattern);
            cellName = pFileBase(idxStart+1:idxEnd);

            % Decide on the actual suffix
            if iPFile <= 2
                % Use the cellName and suffix for the first two p files
                suffixesActual{iPFile} = ['_', cellName, suffixes{iPFile}]; 
            else
                % Use the entire pfile base for other p files
                % Note: For future use
                suffixesActual{iPFile} = ['_', pFileBase];
            end

            % Load p file
            pfiledata = importdata(fullfile(subDirPath, pFileName));    

            % Note: paramNames should be in the fisrt column
            paramNames{iPFile} = pfiledata(:, 1);
            
            % Note: paramValues should be in the second column
            paramValues{iPFile} = cell2mat(pfiledata(:, 2));
        end

        % Compare the parameters and plot output in outFolder
        m3ha_compare_neuronparams2 (paramValues, paramNames, suffixesActual, ...
                                    'OutFolder', outFolder);
    else
        fprintf('pFileNames doesn''t exist for %s!\n', subDirPath);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% clear paramValues paramNames;
% uncomment if running not in parfor
% clear pFileNames; % uncomment if running not in parfor

%OutfolderDefault = '/media/shareX/share/Brian/compare_neuronparams2/';
%OutfolderDefault = '/media/adamX/m3ha/optimizer4gabab/';

addRequired(iP, 'identifier', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

m3ha_compare_neuronparams2 (paramValues, paramNames, suffixes, ...
                            'OutFolder', [outFolder subDirNames{iSubDirs}]);

pFileNames = subFiles(cellfun(@(x) any(regexp(x, '\.p$')), subFiles) == 1);

addParameter(iP, 'Suffixes', suffixesDefault, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));

paramNames = {};
paramValues = {};

pfiledata = importdata(fullfile(subDirNames{iSubDirs}, pFileNames{iSuffix}));    

% paramNames in first column
sub_paramNames = pfiledata(:,1);

% paramValues in other columns
sub_paramValues = cell2mat(pfiledata(:,2:end));

if isempty(paramNames)
    % first time adding paramNames
    paramNames = sub_paramNames;
else
    % add new paramNames to current cell array
    paramNames = {paramNames sub_paramNames};
end
if isempty(paramValues);
    % first time adding paramValues
    paramValues = sub_paramValues;
else
    % Add new paramValues to current cell array
    paramValues = {paramValues sub_paramValues};
end

dirIdentifier = iP.Results.DirIdentifier;
addParameter(iP, 'DirIdentifier', dirIdentifierDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

dirIdentifierDefault = '';

inputDirectoryPattern = ['*', dirIdentifier, '*'];

m3ha_compare_neuronparams2 (paramValues, paramNames, suffixes, ...
                            'OutFolder', outFolder);

pFileNames2 = ...
    subFileNames(cellfun(@(x, y) any(regexp(x, pfilePattern2)), ...
                                 subFileNames));

% Put all matching pFile names together
pFileNames = [pFileNames1, pFileNames2];

if iPFile == 1
    % Use the cellName and suffix for the first p file
    suffixesActual{iPFile} = ['_', cellName, suffixes{1}]; 
elseif iPFile == 2 && numel(suffixes) >= 2
    % Use the suffix only for the second p file
    suffixesActual{iPFile} = suffixes{2};
else
    % Use the entire pfile base for other p files
    suffixesActual{iPFile} = ['_', pFileBase];
end

pFileNames = union_over_cells(pFileNamesAll, 'SetOrder', 'stable');

%}
