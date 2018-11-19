function fileNamesToTakeOut = m3ha_find_files_to_take_out (varargin)
%% Returns all the file names of files to take out from .png files in 'TAKE_OUT_*' folders of a special cases directory
% Usage: fileNamesToTakeOut = m3ha_find_files_to_take_out (varargin)
% Explanation:
%   Returns all the file names of files to take out from casesDir
%   Note: These files are under directories labelled with 'TAKE_OUT_*' 
%           and were the result of voting by blab
% Outputs:
%       fileNamesToTakeOut  - file names of matfiles to take out
%                           specified as a cell array
% Arguments:
%       varargin    - 'CasesDir' - the directory that contains 
%                                   'TAKE_OUT_*' folders with special cases
%                   must be a directory
%                   default == ~/m3ha/data_dclamp/take4/special_cases
%
% Used by:
%       cd/m3ha_select_sweeps_to_fit.m

% File History:
% 2018-11-19 Moved from m3ha_select_sweeps_to_fit.m
% 

%% Hard-coded parameters
takeOutPattern = 'TAKE_OUT_*';
scaledSuffix = '_scaled';
oldSuffix = '_old';
matExtension = '.mat';
pngExtension = '.png';

%% Default values for optional arguments
casesDirDefault = '~/m3ha/data_dclamp/take4/special_cases';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'CasesDir', casesDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, varargin{:});
casesDir = iP.Results.CasesDir;

%% Preparation
% Check if casesDir exists
if ~isfolder(casesDir)
    fprintf('%s does not exist!\n', casesDir);
    fileNamesToTakeOut = {};
    return;
end

% Construct the _scaled.png file ending
scaledEnding = [scaledSuffix, pngExtension];

%% Do the job
% Get all folders with sweeps to take out from fitting
allFoldersToTakeOut = dir(fullfile(casesDir, takeOutPattern));

% Start a counter for files to take out
ct = 0;

% Go through all folders starting with 'TAKE_OUT_*'
for folder = allFoldersToTakeOut'
    % Extract the current folder name
    folderName = folder.name;

    % Get all .png files in this directory
    potentialFilesToTakeOut = ...
        dir(fullfile(casesDir, folderName, ['*', pngExtension]));

    % Go through each .png file
    for file = potentialFilesToTakeOut'
        % Extract the current file name
        fileName = file.name;

        % Count each sweep only once by using only the _scaled file
        if contains(fileName, scaledSuffix) && ...
            ~contains(fileName, oldSuffix)
            % Increment the counter
            ct = ct + 1;

            % Replace the file ending with the matfile ending
            %   Note: this is used for storing sweep info
            matFileName = replace(fileName, scaledEnding, matExtension);

            % Add matfile name for the sweep to take out
            fileNamesToTakeOut{ct} = matFileName;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

matFileName = strrep(fileName, scaledEnding, matExtension);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%