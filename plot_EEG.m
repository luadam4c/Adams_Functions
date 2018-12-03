function plot_EEG(varargin)
%% Plots EEG traces from .abf file(s) or all .abf files in a directory 
% Usage: plot_EEG(varargin)
% Explanation:
%       TODO
% Example(s):
%       plot_EEG('WAGS04_30_2018_cage1.abf');
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       varargin    - 'Directory': the name of the directory containing 
%                                   the .abf files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileNames': names of .abf files to load
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from pwd
%                   - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'OutFolder': the name of the directory in which 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == same as 'Directory'
%                   
% Requires:
%       cd/plot_swd_raster.m
%       cd/plot_traces_EEG.m


% File History:
%   2018-11-26 - Created
% 

%% Default values for optional arguments
directoriesDefault = '';    % set later
fileNamesDefault = {};      % detect from pwd by default
verboseDefault = true;
outFolderDefault = '';      % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoriesDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
verbose = iP.Results.Verbose;
outFolder = iP.Results.OutFolder;

%% Preparation
% Decide on the directory
if isempty(directory)
    % Use the present working directory
    directory = pwd;
end

% Decide on the files to use
if isempty(fileNames)
    % Find all .abf files in the directory
    [~, fileNames] = all_files('Directory', directory, ...
                                'Extension', '.abf', ...
                                'Verbose', verbose);

    % Return usage message if no .abf files found
    if isempty(fileNames)
        fprintf('No .abf files found in %s!\n', directory);
        fprintf('Type ''help %s'' for usage\n', mfilename);
        return
    end
    
    % If there is only one file, change to a character array
    if numel(fileNames) == 1
        fileNames = fileNames{1};
    end
else
    % Display warning if the directory is not empty
    % TODO
end

% Decide on the output folder
if isempty(outFolder)
    outFolder = directory;
end
if verbose
    fprintf('Outfolder is %s ...\n', outFolder);
end

% Check if needed output directory exist
check_dir(outFolder, 'Verbose', verbose);

%% Plot traces for each .abf file
if ischar(fileNames)
    % Plot all traces for this .abf file
    plot_traces_EEG(fileNames, 'Verbose', verbose);
elseif iscell(fileNames)
    % Count the number of files
    nAbfFiles = numel(fileNames);

    % Run through each .abf file
    for iAbfFile = 1:nAbfFiles
        % Get the current .abf file name
        abfFileName = fileNames{iAbfFile};

        % Plot all traces for this .abf file
        plot_traces_EEG(abfFileName, 'Verbose', verbose);
    end
end

%% Plot raster plots
% Plot a raster plot for all SWD files
% TODO: Pass in tables
% TODO: Pass in time limits detected from data files
plot_swd_raster;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       abfFileName - file name of the abf file
%                       could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%