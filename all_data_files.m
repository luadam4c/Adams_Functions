function [dataFiles, dataPaths, extension, message] = all_data_files (varargin)
%% Looks for data files in a directory according to either extensionUser or going through a list of possibleExtensions
% Usage: [dataFiles, dataPaths, extension, message] = all_data_files (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [dataFiles, dataPaths, extension, message] = all_data_files;
%
% Outputs:
%       TODO
%
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ShowFlag': whether to show message
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'WarnFlag': whether to warn if no files found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ExtensionUser'- data file extension provided by user
%                                       or 'auto' for automatic detection
%                   must be a string scalar or a character vector
%                   default == 'auto'
%                   - 'PossibleExtensions': possible data file extensions
%                                           ordered by precedence
%                   must be a string scalar or a character vector
%                   default == {'abf', 'mat', 'txt'}
%                   - 'ExcludedStrings': excluded strings from data file names
%                   must be a cell array of character vectors
%                   default == {}
%                   - Any other parameter-value pair for all_files()
%
% Requires:
%       cd/all_files.m
%       cd/construct_and_check_fullpath.m
%       cd/isemptycell.m
%
% Used by:
%       cd/minEASE.m
%       cd/combine_sweeps.m

% File History:
%   2018-01-29 Moved from cd/minEASE.m
%   2018-01-29 Added the case where extensionUser is not recognized 
%   2018-01-29 Added Keyword as an optional argument
%   2018-02-14 Added ExcludedStrings as an optional argument
%   2020-08-26 Now uses all_files.m and sorts files by datenum
%   2020-08-26 Made all arguments optional
%   2020-08-26 Redefined outputs to be consistent with all_files.m
%   2020-08-26 Removed 'FileIdentifier' and passes extra arguments
%                tod all_files.m instead
%   2020-08-26 Added 'Verbose' and 'ShowFlag' as optional arguments
%   2020-08-27 Added 'WarnFlag' as an optional argument

%% Default values for optional arguments
verboseDefault = false;         % don't print to standard output by default
showFlagDefault = true;         % show message by default
warnFlagDefault = false;        % don't warn if no files found by default
directoryDefault = '';          % construct_and_check_fullpath('') == pwd
extensionUserDefault = 'auto';
possibleExtensionsDefault = {'abf', 'mat', 'txt'};     
excludedStringsDefault = {};    % no excluded strings by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ShowFlag', showFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'WarnFlag', warnFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['Directory must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'ExtensionUser', extensionUserDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PossibleExtensions', possibleExtensionsDefault, ...
    @iscellstr);
addParameter(iP, 'ExcludedStrings', excludedStringsDefault, ...
    @iscellstr);

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
showFlag = iP.Results.ShowFlag;
warnFlag = iP.Results.WarnFlag;
directory = iP.Results.Directory;
extensionUser = iP.Results.ExtensionUser;
possibleExtensions = iP.Results.PossibleExtensions;
excludedStrings = iP.Results.ExcludedStrings;

% Keep unmatched arguments for the all_files() function
otherArguments = iP.Unmatched;

%% Preparation
% Make sure the directory is an existing full path
directory = construct_and_check_fullpath(directory);

%% Find all data files
% Extract the number of possible data types
nDataTypes = numel(possibleExtensions);

% Find data files according to data type
switch extensionUser
case 'auto'                         % if automatically detecting data types
    % Iterate through possibleExtensions to look for possible data type
    for iDataType = 1:nDataTypes
        tempExtension = possibleExtensions{iDataType};
        [dataFiles, dataPaths] = ...
            find_valid_files(directory, tempExtension, ...
                                excludedStrings, warnFlag, otherArguments);
        if numel(dataFiles) > 0
            extension = tempExtension;
            icon = 'none';
            message = {sprintf(['The .%s files in this directory ', ...
                                'will be used as data:'], extension), ...
                        sprintf('%s', directory)};
            break;
        end
    end
    if numel(dataFiles) == 0
        extension = '';
        icon = 'warn';
        message = {'There are no acceptable data files in this directory:', ...
                        sprintf('%s', directory)};
    end
case possibleExtensions              % if user provided a possible data type
    extension = extensionUser;
    [dataFiles, dataPaths] = ...
            find_valid_files(directory, extensionUser, ...
                                excludedStrings, warnFlag, otherArguments);
    if numel(dataFiles) == 0
        icon = 'warn';
        message = {sprintf('There are no .%s files in this directory:', ...
                            extensionUser), sprintf('%s', directory)};
    else
        icon = 'none';
        message = {sprintf(['The .%s files in this directory ', ...
                            'will be used as data:'], extension), ...
                    sprintf('%s', directory)};
    end
otherwise
    extension = '';
    dataFiles = [];
    icon = 'none';

    % Start message
    message = {sprintf('The data type %s is not recognized!', extensionUser), ...
                sprintf('The possible data types are: %s', ...
                        strjoin(possibleExtensions, ', '))};                     
  
    % End message
    message = [message, 'You could also say ''auto''.\n'];
end

%% Print appropriate message
if showFlag
    print_or_show_message(message, 'Icon', icon, ...
                            'MessageMode', 'show', 'Verbose', verbose);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataFiles, dataPaths] = ...
                find_valid_files (directory, extension, ...
                                    excludedStrings, warnFlag, otherArguments)

% Find all files with the given name
[dataFiles, dataPaths] = ...
    all_files('Directory', directory, 'SortBy', 'datenum', ...
                'Extension', extension, 'WarnFlag', warnFlag, otherArguments);

% Exclude invalid entries by testing the date field
toKeep = ~isemptycell({dataFiles.date});
dataFiles = dataFiles(toKeep);
dataPaths = dataPaths(toKeep);

% Exclude entries with an excluded string in the file name
for iString = 1:numel(excludedStrings)
    if ~isempty(dataFiles)       % if dataFiles not already empty
        % Get this excluded string
        string = excludedStrings{iString};
        
        % Get all the data file names as a cell array
        allNames = {dataFiles.name};
        
        % Determine for each file whether you cannot find the string in the name
        doesNotContainString = isemptycell(strfind(allNames, string));
        
        % Restrict to those files 
        dataFiles = dataFiles(doesNotContainString);
        dataPaths = dataPaths(doesNotContainString);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
