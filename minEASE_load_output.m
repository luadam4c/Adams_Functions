function [outputStructArray, outputCellArrays] = ...
                minEASE_load_output (varargin)
%% Load the output and relevant parameters from a minEASE output directory into an output structure
% Usage: [outputStructArray, outputCellArrays] = ...
%               minEASE_load_output (varargin)
% Examples:
%       minEASE_load_output
%       minEASE_load_output('LoadMode', 'combined')
%       minEASE_load_output('LoadMode', 'individual')
%       minEASE_load_output('OutputDir', someDir)
%       minEASE_load_output('OutputFiles', {file1, file2})
%       minEASE_load_output('OutputDir', someDir, 'OutputFiles', {file1, file2})
% Outputs:
%       outputStructArray - a structure array with the following fields:
%                   'outputPath', 'fileIdentifier', 'directionLabel', ...
%                   'sweepNumber', 'timeUnits', 'fileType', ...
%                   'eventInfo', 'eventClass', 'isChecked', ...
%                   'nSamples', 'siMs', 'prevSweepsDuration', ...
%                   'outputLabel'
%       outputCellArray   - a structure with the following fields 
%                           (all cell arrays):
%                   'outputPathAll', 'fileIdentifierAll', 'directionLabelAll', ...
%                   'sweepNumberAll', 'timeUnitsAll', 'fileTypeAll', ...
%                   'eventInfoAll', 'eventClassAll', 'isCheckedAll', ...
%                   'nSamplesAll', 'siMsAll', 'prevSweepsDurationAll', ...
%                   'outputLabelAll'
% Arguments:
%       varargin    - 'OutputFiles': output file names from a minEASE run
%                   must be a valid file or a cell array of valid files
%                       under outputDir
%                   default == either '*_ALL_output*samples*.mat'
%                               or '*_output*.mat' from the output directory
%                   - 'OutputDir': output directory from a minEASE run
%                   must be a valid directory
%                   default == pwd or the directory containing the output files
%                   - 'LoadMode': the type of output file to load from
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'combined'   - Load from '*_ALL_output*samples*.mat'
%                       'individual' - Load from '*_output*.mat'
%                       'auto'       - Load from whatever is available,
%                                       with '*_ALL_output*samples*.mat'
%                                       having priority
%                   default == 'auto'              
%
% Requires:
%       cd/extract_from_minEASE_outputfilename.m
%
% Used by:
%       cd/minEASE_filter_output.m

% File History:
%   2018-07-26 Modified from /home/Matlab/Marks_Functions/paula/Oct2017/loadAllTheData.m
%   2018-08-03 Updated usage of minEASE_extract_from_output_filename.m
%   2018-08-03 Changed sweepLabel -> outputLabel
%   2018-08-04 Added 'OutputFiles' as an optional argument
%               and allowed loading of files from different directories
%   2018-08-05 Now outputs both a struct array and a struct of cell arrays
%   2018-09-21 Switched back to '*_output*.mat' for individualOutputPattern
% TODO: Load from csv files if matfiles are not available
% 

%% Hard-coded constants
% To be consistent with minEASE.m
individualOutputPattern = '*_output*.mat';
combinedOutputPattern = '*_ALL_output*samples*.mat';

% Valid load modes
validLoadModes = {'combined', 'individual', 'auto'};

% Fields of the output structure
outputFields = {'outputPath', 'fileIdentifier', 'directionLabel', ...
                'sweepNumber', 'timeUnits', 'fileType', ...
                'eventInfo', 'eventClass', 'isChecked', ...
                'nSamples', 'siMs', 'prevSweepsDuration', 'outputLabel'};

%% Default values for optional arguments
outputFileNamesDefault = {};            % to be set later
outputDirDefault = '';                  % to be set later
loadModeDefault = 'auto';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputFiles', outputFileNamesDefault, ...
    @(x) assert(ischar(x) || iscell(x) && (all(cellfun(@ischar, x)) ...
        || all(cellfun(@isstring, x))) || isstring(x) , ...
        ['OutputFileNames must be either a string/character array ', ...
            'or a cell array of strings/character arrays!']));
addParameter(iP, 'OutputDir', outputDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'LoadMode', loadModeDefault, ...
    @(x) any(validatestring(x, validLoadModes)));

% Read from the Input Parser
parse(iP, varargin{:});
outputFileNames = iP.Results.OutputFiles;
outputDir = iP.Results.OutputDir;
loadMode = validatestring(iP.Results.LoadMode, validLoadModes);

%% Decide on the output files to load
% Initialize an empty structure in case function exits early
outputStructArray = struct;
outputCellArrays = struct;

% Check if user-defined output directory is valid
if ~isempty(outputDir)
    % The directory must exist
    if ~isfolder(outputDir)
        fprintf('%s does not exist or is not readable!\n', outputDir);
        return;
    end
end

% If the output file names provided is not a cell array, convert to a cell array
if ~iscell(outputFileNames)
    outputFileNames = {outputFileNames};
end

% Make sure the cell array is a column array
outputFileNames = outputFileNames(:);

% If user defined output file names, use them. 
%   Otherwise choose the output files base on load mode
if ~isempty(outputFileNames)
    % Break the file names apart
    [outputDirs, outputFileBases, outputFileExts] = ...
        cellfun(@fileparts, outputFileNames, 'UniformOutput', false);

    % Get unique entries in outputDirs
    uniqueOutputDirs = unique(outputDirs);

    % Decide on outputDir and outputPathAll
    % Two cases: more than one or exactly one unique entry
    if numel(uniqueOutputDirs) > 1    
        % If there is more than one entry in outputDirs,
        %   check if outputDir is empty
        if ~isempty(outputDir)
            fprintf('At least one output file does not exist in %s!\n', ...
                    outputDir);
            return;
        else
            % Keep the outputDir empty
        end

        % Use outputFileNames as the full paths to the output files
        outputPathAll = outputFileNames;
    else
        % If there is only one entry in outputDirs,
        %   call that outputDirsCommon (could be empty)
        outputDirsCommon = uniqueOutputDirs{1};

        % Make sure outputDir matches outputDirsCommon (could be empty)
        if ~isempty(outputDir) 
            if ~strcmp(outputDir, outputDirsCommon)
                fprintf(['The directories of output files is %s, ', ...
                        'but the given output directory is %s!\n'], ...
                        outputDirsCommon, outputDir);
                return;
            else
                % The output is already outputDirsCommon
            end
        else
            % Make outputDirsCommon (could still be empty) the output directory
            outputDir = outputDirsCommon;
        end

        % If outputDir is empty, use the present working directory
        if isempty(outputDir)
            outputDir = pwd;
        end

        % Construct the full paths to the output files
        outputPathAll = cellfun(@(x, y) fullfile(outputDir, [x, y]), ...
                                outputFileBases, outputFileExts, ...
                                'UniformOutput', false);
    end

    % Check whether each element is a valid file
    isFile = cellfun(@isfile, outputPathAll);

    % If not all of them are valid files, print message
    if ~all(isFile)
        for iFile = find(~isFile)
            fprintf('%s does not exist or is not readable!\n', ...
                    outputPathAll{iFile});
        end
        return
    end
else
    % If outputDir is empty, use the present working directory
    if isempty(outputDir)
        outputDir = pwd;
    end

    % Choose the output files based on the load mode
    outputFileNames = choose_output_files (outputDir, loadMode, ...
                        combinedOutputPattern, individualOutputPattern);

    % Make sure the cell array is a column array
    outputFileNames = outputFileNames(:);
                    
    % Construct full paths
    outputPathAll = cellfun(@(x) fullfile(outputDir, x), ...
                        outputFileNames, 'UniformOutput', false);
end

% Count the number of files
nFiles = numel(outputPathAll);

% See if the load mode agrees with the number of files
if strcmp(loadMode, 'combined')
    if nFiles > 1
        fprintf('Cannot load more then one file in ''combined'' mode!\n');
        fprintf('Please run minEASE_combine_events.m!\n');
        return
    end
end

%% Load output matfiles in the directory
% Load all matfiles
fileIdentifierAll = cell(nFiles, 1);
directionLabelAll = cell(nFiles, 1);
sweepNumberAll = cell(nFiles, 1);
timeUnitsAll = cell(nFiles, 1);
fileTypeAll = cell(nFiles, 1);
eventInfoAll = cell(nFiles, 1);
eventClassAll = cell(nFiles, 1);
isCheckedAll = cell(nFiles, 1);
nSamplesAll = cell(nFiles, 1);
siMsAll = cell(nFiles, 1);
prevSweepsDurationAll = cell(nFiles, 1);
outputLabelAll = cell(nFiles, 1);
for iFile = 1:nFiles
    % Get the name of the output matfile
    %   Note: this should be in the format:
    %   [fileIdentifier, '_', directionLabel, ...
    %       '_Swp', num2str(iSwp)]_output.mat
    thisFilePath = outputPathAll{iFile};

    % Load info from both the file name and the file
    [fileIdentifierAll{iFile}, directionLabelAll{iFile}, ...
        sweepNumberAll{iFile}, timeUnitsAll{iFile}, fileTypeAll{iFile}, ...
        eventInfoAll{iFile}, eventClassAll{iFile}, isCheckedAll{iFile}, ...
        nSamplesAll{iFile}, siMsAll{iFile}, ...
        prevSweepsDurationAll{iFile}, outputLabelAll{iFile}] = ...
        minEASE_load_output_helper(thisFilePath);
end

% Place variables in the outputCellArrays structure
outputCellArrays.outputPathAll = outputPathAll;
outputCellArrays.fileIdentifierAll = fileIdentifierAll;
outputCellArrays.directionLabelAll = directionLabelAll;
outputCellArrays.sweepNumberAll = sweepNumberAll;
outputCellArrays.timeUnitsAll = timeUnitsAll;
outputCellArrays.fileTypeAll = fileTypeAll;
outputCellArrays.eventInfoAll = eventInfoAll;
outputCellArrays.eventClassAll = eventClassAll;
outputCellArrays.isCheckedAll = isCheckedAll;
outputCellArrays.nSamplesAll = nSamplesAll;
outputCellArrays.siMsAll = siMsAll;
outputCellArrays.prevSweepsDurationAll = prevSweepsDurationAll;
outputCellArrays.outputLabelAll = outputLabelAll;

% Place everything in an overall cell array
masterCellArray = [outputPathAll, fileIdentifierAll, directionLabelAll, ...
                    sweepNumberAll, timeUnitsAll, fileTypeAll, ...
                    eventInfoAll, eventClassAll, isCheckedAll, ...
                    nSamplesAll, siMsAll, prevSweepsDurationAll, ...
                    outputLabelAll];

% Convert everything to a structure array
outputStructArray = cell2struct(masterCellArray, outputFields, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputFileNames = choose_output_files (outputDir, loadMode, ...
                            combinedOutputPattern, individualOutputPattern)

switch loadMode
case {'auto', 'combined'}
    % Load from the combined output file if available
    allFiles = dir(fullfile(outputDir, combinedOutputPattern));

    % Make sure there is exactly one such file
    nAllFiles = numel(allFiles);
    if nAllFiles > 1
        fprintf(['There should be only one file with ', ...
                    'the pattern ''%s'' in the directory %s!\n'], ...
                combinedOutputPattern, outputDir);
        return;
    end

    % If there is no such file, see if there are 
    %   more than one output matfiles
    if nAllFiles < 1
        % Find all output matfiles in the directory
        outputFiles = dir(fullfile(outputDir, individualOutputPattern));

        % Get the output file names
        outputFileNames = {outputFiles.name};

        % Count the number of output files
        nOutputFiles = numel(outputFileNames);

        % Decide based on the number of output files
        if nOutputFiles > 1
            % If there is more than one file, decide based on load mode
            if strcmp(loadMode, 'combined')
                % Ask the user to combine the output files
                fprintf(['There is more than one output file ', ...
                            'with the pattern ''%s'' in the directory %s!\n'], ...
                        individualOutputPattern, outputDir);
                fprintf('Please run minEASE_combine_events.m!\n');
                return
            else
                % Use those files
            end
        elseif nOutputFiles < 1
            % If there is no output file, return message
            fprintf(['There is no output file ', ...
                        'with the pattern ''%s'' in the directory %s!\n'], ...
                    individualOutputPattern, outputDir);
            return
        else
            % If there is exactly one output file,
            %   use it as the matfile to load
        end
    else
        % Use the combined output file as the matfile to load
        outputFileNames = {allFiles.name};
    end
case 'individual'
    % Find all output matfiles in the directory
    outputFiles = dir(fullfile(outputDir, individualOutputPattern));

    % Get all output file names
    outputFileNames = {outputFiles.name};
otherwise 
    error('Load mode ''%s'' unrecognized!', loadMode);       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fileIdentifier, directionLabel, sweepNumber, timeUnits, fileType, eventInfo, eventClass, isChecked, nSamples, siMs, prevSweepsDuration, outputLabel] = minEASE_load_output_helper(outputPath)
%% Loads info from both the file name and the file

% Extract the output file name
[~, outputFileBase, outputFileExt] = fileparts(outputPath);
outputFileName = [outputFileBase, outputFileExt];

% Extract the fileIdentifier, directionLabel, sweepNumber,
%   timeUnits, fileType
infoStruct = minEASE_extract_from_output_filename(outputFileName);
fileIdentifier = infoStruct.fileIdentifier;
directionLabel = infoStruct.directionLabel;
sweepNumber = infoStruct.sweepNumber;
timeUnits = infoStruct.timeUnits;
fileType = infoStruct.fileType;

% Read in the v7.3 output matfile
output = load(fullfile(outputPath));

% Extract the variables
if isfield(output, 'eventInfo')  
    eventInfo = output.eventInfo;
else
    eventInfo = NaN;
end
if isfield(output, 'eventClass')  
    eventClass = output.eventClass;
else
    eventClass = NaN;
end
if isfield(output, 'isChecked')  
    isChecked = output.isChecked;
else
    isChecked = NaN;
end
if isfield(output, 'nSamples')  
    nSamples = output.nSamples;
else
    nSamples = NaN;
end
if isfield(output, 'siMs')  
    siMs = output.siMs;
else
    siMs = NaN;
end
if isfield(output, 'prevSweepsDuration')  
    prevSweepsDuration = output.prevSweepsDuration;
else
    prevSweepsDuration = NaN;
end
if isfield(output, 'outputLabel')  
    outputLabel = output.outputLabel;
elseif isfield(output, 'sweepLabel')  
    outputLabel = output.sweepLabel;
else
    outputLabel = '';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Read in the params matfile
params = load(fullfile(outputDir, paramsFileName));

if isfield(params, 'nSamples')  
    nSamples = params.nSamples;
else
    nSamples = [];
end
if isfield(params, 'siMs')  
    siMs = params.siMs;
else
    siMs = [];
end
if isfield(params, 'prevSweepsDuration')  
    prevSweepsDuration = params.prevSweepsDuration;
else
    prevSweepsDuration = [];
end

individualParamsPattern = '*_params.mat';
combinedParamsPattern = '*_ALL_params.mat'; % TODO
% Get the name of the parameters matfile
paramsFileName = strrep(outputFileName, individualOutputPattern(2:end), ...
                        individualParamsPattern(2:end));
% Get the name of the parameters matfile
paramsFileName = strrep(outputFileName, combinedOutputPattern(2:end), ...
                        combinedParamsPattern(2:end));
% Get the name of the parameters matfile
paramsFileName = strrep(outputFileName, individualOutputPattern(2:end), ...
                        individualParamsPattern(2:end));
output.paramsFileName = paramsFileName;
paramsFileNames = cell(nFiles, 1);
output.paramsFileNames = paramsFileNames;

individualOutputPattern = '*_output.mat';
combinedOutputPattern = '*_ALL_output_samples.mat';

% Get full path to the output file
outputFile = fullfile(outputDir, outputFileName);

% Deal with a single file and many files differently
if ischar(outputFile) && ~isfile(outputFile)
    fprintf('%s does not exist or is not readable!\n', outputFile);
    return
elseif iscell(outputFile)
end

fprintf('The chosen output file is %s\n', outputFile);

if strcmp(loadMode, 'combined') || ...
    (strcmp(loadMode, 'auto') && nFiles == 1)
    % If there is more than one file provided,
    %   ask user to combine them
    if nFiles > 1
        fprintf('Cannot load more then one file in ''combined'' mode!\n');
        fprintf('Please run minEASE_combine_events.m!\n');
        return
    end

    % Otherwise, load that file
    outputPath = outputPaths{1};
    [fileIdentifier, directionLabel, ...
        sweepNumber, timeUnits, fileType, ...
        eventInfo, eventClass, isChecked, ...
        nSamples, siMs, prevSweepsDuration, outputLabel] = ...
        minEASE_load_output_helper(outputPath);

    % Place variables in the output structure
    output.outputPath = outputPath;
    output.fileIdentifier = fileIdentifier;
    output.directionLabel = directionLabel;
    output.sweepNumber = sweepNumber;
    output.timeUnits = timeUnits;
    output.fileType = fileType;
    output.eventInfo = eventInfo;
    output.eventClass = eventClass;
    output.isChecked = isChecked;
    output.nSamples = nSamples;
    output.siMs = siMs;
    output.prevSweepsDuration = prevSweepsDuration;
    output.outputLabel = outputLabel;

elseif strcmp(loadMode, 'individual') || ...
        (strcmp(loadMode, 'auto') && nFiles > 1)
    % Load all matfiles
    fileIdentifiers = cell(nFiles, 1);
    directionLabels = cell(nFiles, 1);
    sweepNumbers = cell(nFiles, 1);
    timeUnitsAll = cell(nFiles, 1);
    fileTypes = cell(nFiles, 1);
    eventInfos = cell(nFiles, 1);
    eventClasses = cell(nFiles, 1);
    isCheckeds = cell(nFiles, 1);
    nSamplesAll = cell(nFiles, 1);
    siMsAll = cell(nFiles, 1);
    prevSweepsDurationAll = cell(nFiles, 1);
    outputLabelAll = cell(nFiles, 1);
    for iFile = 1:nFiles
        % Get the name of the output matfile
        %   Note: this should be in the format:
        %   [fileIdentifier, '_', directionLabel, ...
        %       '_Swp', num2str(iSwp)]_output.mat
        thisFilePath = outputPaths{iFile};

        % Load info from both the file name and the file
        [fileIdentifiers{iFile}, directionLabels{iFile}, ...
            sweepNumbers{iFile}, timeUnitsAll{iFile}, fileTypes{iFile}, ...
            eventInfos{iFile}, eventClasses{iFile}, isCheckeds{iFile}, ...
            nSamplesAll{iFile}, siMsAll{iFile}, ...
            prevSweepsDurationAll{iFile}, outputLabelAll{iFile}] = ...
            minEASE_load_output_helper(thisFilePath);
    end

    % Place variables in the output structure
    output.outputPaths = outputPaths;
    output.fileIdentifiers = fileIdentifiers;
    output.directionLabels = directionLabels;
    output.sweepNumbers = sweepNumbers;
    output.timeUnitsAll = timeUnitsAll;
    output.fileTypes = fileTypes;
    output.eventInfos = eventInfos;
    output.eventClasses = eventClasses;
    output.isCheckeds = isCheckeds;
    output.nSamplesAll = nSamplesAll;
    output.siMsAll = siMsAll;
    output.prevSweepsDurationAll = prevSweepsDurationAll;
    output.outputLabelAll = outputLabelAll;
else
    error('The load mode %s is not recognized!\n', loadMode);
end

% individualOutputPattern = '*_output.mat';

%}
