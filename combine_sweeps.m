function [allData, timeVec] = combine_sweeps (varargin)
%% Combines sweeps that begin with expLabel in dataDirectory under dataMode
% Usage: [allData, timeVec] = combine_sweeps (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       TODO
%
% Arguments:    
%       varargin    - 'DataDirectory'- dataDirectory containing files to combine
%                   must be a string scalar or a character vector
%                   default == set in all_data_files.m
%                   - 'FilePaths': names of '.mat' output files to combine
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from dataDirectory
%                   - 'ExpLabel'    - experiment label for files to combine
%                   must be a string scalar or a character vector
%                   default == '' 
%                   - 'FileIdentifier': data file identifier (may be empty)
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'OutputLabel' - experiment label for output file names
%                   must be a string scalar or a character vector
%                   default == expLabel if provided and 
%                                   dataDirectory name otherwise
%                   - 'SweepNumbers' - the sweep numbers to combine
%                   must be a numeric vector or 'all'
%                   default == 'all'
%                   - 'Extension': input data type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'abf'   - AXON binary files
%                       'mat'   - MATLAB data files
%                       'txt'   - text files
%                       'auto'  - automatically detect from dataDirectory
%                   default == 'auto'
%                   - 'DataMode': input data mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '2d'    - a single sweep per file
%                       '3d'    - multiple sweeps per file
%                       'auto'  - automatically detect from dataDirectory
%                   default == 'auto'
%                   - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%                   - 'Verbose' - whether to print to standard output
%                                   regardless of message mode
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/all_data_files.m
%       cd/construct_fullpath.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/force_row_cell.m
%       cd/identify_channels.m
%       cd/locate_functionsdir.m
%       cd/print_or_show_message.m
%       /home/Matlab/Downloaded_Functions/abf2load.m
%
% Used by:    
%       cd/minEASE.m
%       /home/Matlab/Katies_Functions/Glucose/cell_attached_minEASE_filtered.m
%       /home/Matlab/Katies_Functions/Glucose/GABAB_puff.m
%       /home/Matlab/Katies_Functions/Glucose/GABAB_puff_test.m
%       /home/Matlab/Katies_Functions/Glucose/loadcell_attached_TimeSeries.m

% File History:
% 2017-07-25 Created by AL
% 2017-10-15 Added success message
% 2018-01-24 Added isdeployed
% 2018-01-29 Now uses all_data_files.m
% 2018-01-29 Added extension as an optional parameter-value pair argument
% 2018-02-02 Added showMessage as an optional parameter-value pair argument
% 2018-02-02 Now uses print_or_show_message.m for output
% 2018-02-07 MD - Changed usage of print_or_show_message()
% 2018-02-27 AL - Changed showMessages to messageMode with possible values:
%                   'wait', 'show', 'none'
% 2018-03-02 MD - Defined verbose parameter for print_or_show_message
% 2018-08-12 AL - Made 'DataDirectory', 'ExpLabel', 'DataMode', etc.
%                   optional arguments
% 2018-08-12 AL - Now detects data mode automatically
% 2018-09-23 AL - Added FileIdentifier as an optional argument
% 2020-08-27 AL - Added 'FilePaths' as an optional argument
% 

%% Hard-coded constants
US_PER_MS = 1000;

%% Hard-coded parameters
possibleExtensions = {'abf', 'mat', 'txt'};     
                    % Precedence: .abf > .mat > .txt
possibleDataModes = {'2d', '3d'};
validMessageModes = {'wait', 'show', 'none'};
expMode = 'patch';              % TODO: Make an optional argument

%% Default values for optional arguments
dataDirectoryDefault = pwd;     % use the present working directory by default
filePathsDefault = {};          % set later
expLabelDefault = '';           % disregard any experiment label by default
fileIdentifierDefault = '';     % no file identifier by default
outputLabelDefault = '';        % (will be changed later)
sweepNumbersDefault = 'all';    % combine all sweeps by default
extensionDefault     = 'auto';  % to detect input data type 
                                %   from possibleExtensions
dataModeDefault     = 'auto';   % to detect input data type
messageModeDefault  = 'none';   % print to standard output by default
verboseDefault = false;         % default: Program does not print message
                                %   even if message box is shown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('abf2load.m', 'file') ~= 2 && ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for abf2load.m, abfload.m
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions'));
end
                                                
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataDirectory', dataDirectoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FilePaths', filePathsDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ExpLabel', expLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileIdentifier', fileIdentifierDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutputLabel', outputLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SweepNumbers', sweepNumbersDefault);
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) any(validatestring(x, possibleExtensions)));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) any(validatestring(x, possibleDataModes)));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
dataDirectory = iP.Results.DataDirectory;
filePaths = iP.Results.FilePaths;
expLabel = iP.Results.ExpLabel;
fileIdentifier = iP.Results.FileIdentifier;
outputLabel = iP.Results.OutputLabel;
sweepNumbers = iP.Results.SweepNumbers;
dataModeUser = validatestring(iP.Results.Extension, ...
                    [possibleDataModes, {'auto'}]);
extensionUser = validatestring(iP.Results.Extension, ...
                    [possibleExtensions, {'auto'}]);
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);
verbose = iP.Results.Verbose;

% Extract from experiment label if requested
if strcmpi(fileIdentifier, 'fromExpLabel')
    fileIdentifier = replace(expLabel, {'_IPSC', '_EPSC'}, {'', ''});
end

% Decide on files to combine
if isempty(filePaths)
    % Set default directory
    if isempty(directory)
        directory = pwd;
    end

    % Determine data type, list all .abf, .mat or .txt files from data directory
    [~, filePaths, extension, message] = ...
        all_data_files('Directory', dataDirectory, ...
                        'ExtensionUser', extensionUser, ...
                        'PossibleExtensions', possibleExtensions, ...
                        'Keyword', fileIdentifier, 'ShowFlag', false);
else
    filePaths = construct_fullpath(filePaths, 'Directory', dataDirectory);
    dataDirectory = extract_common_directory(filePaths);
    extensions = extract_fileparts(filePaths, 'Extension');
    extension = extractAfter(extensions{1}, '.');    
    message = [{'Sweeps from: '}, force_row_cell(filePaths)];
end

% Count the number of data files
nDataFiles = numel(filePaths);

% Check if files are found
if nDataFiles <= 0
    if iscell(message)
        message = [message, {'Failed to combine sweeps!'}];
    else
        message = {message, 'Failed to combine sweeps!'};
    end
    mTitle = 'Data files not found';
    % TODO: load custom icon
    icon = 'error';
    print_or_show_message(message, 'MessageMode', messageMode, ...
                            'MTitle', mTitle, 'Icon', icon, ...
                            'Verbose', verbose);
    return;
else
    % TODO: implement icon
    mTitle = 'Data type used';
    % TODO: load custom icon
    icon = 'none';
    print_or_show_message(message, 'MessageMode', messageMode, ...
                            'MTitle', mTitle, 'Icon', icon, ...
                            'Verbose', verbose);
end

% Decide on the data mode and the sampling interval if any
if strcmpi(dataModeUser, 'auto')
    if strcmpi(extension, 'abf')
        % Get the first data file name
        firstDataPath = filePaths{1};

        % Load it
        % TODO: Use parse_abf.m
        [data, siUs, ~] = abf2load(firstDataPath);
%         [data, siUs] = abfload(firstDataPath);

        % Determine the dimensions of the data structure
        nDims = length(size(data));

        % Decide based on the dimensions
        if nDims == 2
            dataMode = '2d';
        elseif nDims == 3
            dataMode = '3d';
        else
            error('data mode unrecognized!');
        end

        % Get the sampling interval in milliseconds
        siMs = siUs / US_PER_MS;
    else
        % The other file types should always be one sweep per file
        dataMode = '2d';

        % TODO: let the user provide a sampling interval
        siMs = [];
    end
else
    % Let the user override the data mode detection
    dataMode = dataModeUser;

    % TODO: let the user provide a sampling interval
    siMs = [];
end

% Concatenate sweep data
allData = [];                       % Initialize combined sweep data
if strcmpi(dataMode, '2d')       % if there is only one sweep per file
    % Set default sweep numbers for one sweep per file
    if ~isnumeric(sweepNumbers) && strcmpi(sweepNumbers, 'all')
        sweepNumbers = 1:nDataFiles;
    end

    % Check if number of abf files is correct
    if nDataFiles < max(sweepNumbers)
        message = 'Too few data files! Failed to combine sweeps!';
        mTitle = 'Combine sweep error';
        icon = 'error';
        print_or_show_message(message, 'MessageMode', messageMode, ...
                                'MTitle', mTitle, 'Icon', icon, ...
                                'Verbose', verbose);
        return;
    end

    % Do for each sweep/file
    for iFile = sweepNumbers
        % Get current data file name
        filePathThis = filePaths{iFile};
        
        % Import data based on data type
        switch extension
        case 'abf'
            % Load abf data for this sweep
            [data, ~] = abf2load(filePathThis);
        case {'mat', 'txt'}
            % Import data (the current vector must be one of the columns)
            %   Only one cell per file!
            data = importdata(filePathThis);            
        otherwise
            message = {sprintf('The data type .%s is not supported yet!', ...
                        extension), 'Failed to combine sweeps!'};
            mTitle = 'Combine sweep error';
            icon = 'error';
            print_or_show_message(message, 'MessageMode', messageMode, ...
                                    'MTitle', mTitle, 'Icon', icon, ...
                                    'Verbose', verbose);
            return;            
        end

        % Identify the current channel from data
        %   To avoid confusion, place only one channel in data
        channelTypes = identify_channels(data, 'ExpMode', expMode);
        if isempty(channelTypes{1})
            idxCurrent = 1;
        else
            idxCurrent = strcmpi('Current', channelTypes);
        end

        % Get current vector (usually pA) for this sweep
        % Note: assume the 2nd dimension is the channel
        current = data(:, idxCurrent);

        % Add to allData
        allData = [allData; current];
    end
elseif strcmpi(dataMode, '3d')   % if there are multiple sweeps per file
    % Check if multiple data files exist
    if nDataFiles > 1
        message = 'Too many data files! Failed to combine sweeps!';
        mTitle = 'Combine sweep error';
        icon = 'error';
        print_or_show_message(message, 'MessageMode', messageMode, ...
                                'MTitle', mTitle, 'Icon', icon, ...
                                'Verbose', verbose);
        return;
    end

    % Get current data file name
    filePathThis = filePaths{1};

    % Load abf data for all sweeps
    [data, ~] = abf2load(abfFileName);
    
    % Find the number of sweeps
    nSwps = size(data, 3);
    
    % Set default sweep numbers for one sweep per file
    if ~isnumeric(sweepNumbers) && strcmpi(sweepNumbers, 'all')
        sweepNumbers = 1:nSwps;
    end

    % Identify the channels in data
    channelTypes = identify_channels(data, 'ExpMode', expMode);
    idxCurrent = strcmpi('Current', channelTypes);
  
    % Do for each sweep
    for iSwp = sweepNumbers
        % Note: assume the 2nd dimension is the channel
        %       and the 3rd dimension is the sweep
        current = data(:, idxCurrent, iSwp);

        % Add to allData
        allData = [allData; current];
    end
end

% Construct a time vector if possible
if ~isempty(siMs)
    % Count the number of samples
    nSamples = length(allData);

    % Construct a time vector in milliseconds
    timeVec = (1:nSamples)' * siMs;
else
    timeVec = [];
end

%% Save output TODO
% Decide on the output label
[~, dataDirectoryName, ~] = fileparts(dataDirectory);
if isempty(outputLabel)
    if ~isempty(expLabel)
        outputLabel = expLabel;
    else
        outputLabel = dataDirectoryName;
    end
end

%% Show success message
message = 'Sweeps successfully combined!!';
mTitle = 'Combine sweep success';
% TODO: load custom icon
icon = 'none';
print_or_show_message(message, 'MessageMode', messageMode, ...
                        'MTitle', mTitle, 'Icon', icon, ...
                        'Verbose', verbose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
