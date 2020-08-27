function [allData, timeVec] = combine_sweeps (varargin)
%% Combines sweeps that begin with expLabel in dataDirectory under dataMode
% Usage: [allData, timeVec] = combine_sweeps (varargin)
% Explanation:
%       TODO
% Outputs:
%       TODO
% Arguments:    
%       varargin    - 'DataDirectory'- dataDirectory containing files to combine
%                   must be a string scalar or a character vector
%                   default == pwd
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
%                   - 'DataType': input data type
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
%       cd/identify_channels.m
%       cd/all_data_files.m
%       cd/locate_functionsdir.m
%       cd/print_or_show_message.m
%       /home/Matlab/Downloaded_Functions/abf2load.m
%
% Used by:    
%       cd/minEASE.m
%       /home/Matlab/Katies_Functions/loadcell_attached_TimeSeries.m
%       /home/Matlab/Katies_Functions/cell_attached_minEASE_filtered.m
%
% File History:
% 2017-07-25 Created by AL
% 2017-10-15 Added success message
% 2018-01-24 Added isdeployed
% 2018-01-29 Now uses all_data_files.m
% 2018-01-29 Added dataType as an optional parameter-value pair argument
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
% 

%% Hard-coded constants
US_PER_MS = 1000;

%% Hard-coded parameters
possibleDataTypes = {'abf', 'mat', 'txt'};     
                    % Precedence: .abf > .mat > .txt
possibleDataModes = {'2d', '3d'};
validMessageModes = {'wait', 'show', 'none'};
expMode = 'patch';              % TODO: Make an optional argument

%% Default values for optional arguments
dataDirectoryDefault = pwd;            % use the present working directory by default
expLabelDefault = '';           % disregard any experiment label by default
fileIdentifierDefault = '';     % no file identifier by default
outputLabelDefault = '';        % (will be changed later)
sweepNumbersDefault = 'all';    % combine all sweeps by default
dataTypeDefault     = 'auto';   % to detect input data type 
                                %   from possibleDataTypes
dataModeDefault     = 'auto';   % to detect input data type
messageModeDefault  = 'none';   % print to standard output by default
verboseDefault = false;         % default: Program does not print message
                                %   even if message box is shown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if ~isdeployed
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
addParameter(iP, 'ExpLabel', expLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileIdentifier', fileIdentifierDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutputLabel', outputLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SweepNumbers', sweepNumbersDefault);
addParameter(iP, 'DataType', dataTypeDefault, ...
    @(x) any(validatestring(x, possibleDataTypes)));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) any(validatestring(x, possibleDataModes)));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
dataDirectory = iP.Results.DataDirectory;
expLabel = iP.Results.ExpLabel;
fileIdentifier = iP.Results.FileIdentifier;
outputLabel = iP.Results.OutputLabel;
sweepNumbers = iP.Results.SweepNumbers;
dataModeUser = validatestring(iP.Results.DataType, ...
                    [possibleDataModes, {'auto'}]);
dataTypeUser = validatestring(iP.Results.DataType, ...
                    [possibleDataTypes, {'auto'}]);
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);
verbose = iP.Results.Verbose;

% Extract from experiment label if requested
if strcmpi(fileIdentifier, 'fromExpLabel')
    fileIdentifier = strrep(strrep(expLabel, '_IPSC', ''), '_EPSC', '');
end

% Determine data type, list all .abf, .mat or .txt files from 
%   data subdirectory
[dataType, allDataFiles, nDataFiles, message] = ...
    all_data_files(dataTypeUser, dataDirectory, possibleDataTypes, ...
                        'FileIdentifier', fileIdentifier);

% Check if files are found
if nDataFiles <= 0
    if iscell(message)
        message = [message, 'Failed to combine sweeps!'];
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
    if strcmpi(dataType, 'abf')
        % Get the first data file name
        firstDataFileName = fullfile(dataDirectory, allDataFiles(1).name);

        % Load it
        [data, siUs, ~] = abf2load(firstDataFileName);
%         [data, siUs] = abfload(firstDataFileName);

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
        dataFileName = fullfile(dataDirectory, allDataFiles(iFile).name);
        
        % Import data based on data type
        switch dataType
        case 'abf'
            % Load abf data for this sweep
            [data, ~] = abf2load(dataFileName);
        case {'mat', 'txt'}
            % Import data (the current vector must be one of the columns)
            %   Only one cell per file!
            data = importdata(dataFileName);            
        otherwise
            message = {sprintf('The data type .%s is not supported yet!', ...
                        dataType), 'Failed to combine sweeps!'};
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
    dataFileName = fullfile(dataDirectory, allDataFiles(1).name);

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

% Find all the .abf files in the dataDirectory with this fileIdentifier
abfFiles = dir(fullfile(dataDirectory, [fileIdentifier, '*.abf']));
nAbfFiles = length(abfFiles);

% Check if files are found
if nAbfFiles <= 0 
    fprintf('Cannot find abf files! Failed to combine sweeps!\n\n');
    return;
end

fprintf('Cannot find data files! Failed to combine sweeps!\n\n');
fprintf('Too few data files! Failed to combine sweeps!');
fprintf(['The data type .%s is not supported yet!\n', ...
        ' Failed to combine sweeps!\n\n'], dataType);
fprintf('Too many data files! Failed to combine sweeps!\n\n');
fprintf('Sweeps successfully combined!!\n\n');

%                   - 'ShowMessage': whether to show messages in messages boxes
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
showMessageDefault  = false;        % print to standard output by default
addParameter(iP, 'ShowMessage', showMessageDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
showMessage = iP.Results.ShowMessage;
        if showMessage
            print_or_show_message(message, 'MessageMode', 'show', ...
                                    'MTitle', mTitle, 'Icon', icon);
        else
            print_or_show_message(message, 'MessageMode', 'none', ...
                                    'MTitle', mTitle, 'Icon', icon);
        end

addpath(fullfile(functionsDirectory, '/Brians_Functions/'));
                                            % for identify_channels.m
addpath(fullfile(functionsDirectory, '/Adams_Functions/'));
                                            % for all_data_files.m

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
