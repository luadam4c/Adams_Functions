function allData = combine_sweeps(dataDirectory, expLabel, dataMode, varargin)
%% Combines sweeps that begin with expLabel in dataDirectory under dataMode
% Usage: allData = combine_sweeps(dataDirectory, expLabel, dataMode, varargin)
% Explanation:
%       TODO
% Outputs:
%       TODO
% Arguments:    
%       TODO
%       varargin    - 'DataType': input data type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'abf'   - AXON binary files
%                       'mat'   - MATLAB data files
%                       'txt'   - text files
%                   default == automatically detect from dataDirectory
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
%       /home/Matlab/Downloaded_Functions/abf2load.m
%       /home/Matlab/Brians_Functions/identify_channels.m
%       /home/Matlab/Adams_Functions/find_data_files.m
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%
% Used by:    
%       /home/Matlab/minEASE/minEASE.m
%
% File History:
% 2017-07-25 Created by AL
% 2017-10-15 Added success message
% 2018-01-24 Added isdeployed
% 2018-01-29 Now uses find_data_files.m
% 2018-01-29 Added dataType as an optional parameter-value pair argument
% 2018-02-02 Added showMessage as an optional parameter-value pair argument
% 2018-02-02 Now uses print_or_show_message.m for output
% 2018-02-07 MD - Changed usage of print_or_show_message()
% 2018-02-27 AL - Changed showMessages to messageMode with possible values:
%                   'wait', 'show', 'none'
% 2018-03-02 MD - Defined verbose parameter for print_or_show_message
% 

%% Hard-coded parameters
possibleDataTypes = {'abf', 'mat', 'txt'};     
                    % Precedence: .abf > .mat > .txt
validMessageModes = {'wait', 'show', 'none'};

%% Default values for optional arguments
sweepNumbersDefault = 'all';
dataTypeDefault     = 'auto';       % to detect input data type 
                                    %   from possibleDataTypes
messageModeDefault  = 'none';       % print to standard output by default
verboseDefault = false;             % default: Program does not print message
                                    %   even if message box is shown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if ~isdeployed
    if exist('/home/Matlab/', 'dir') == 7
        functionsDirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsDirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsDirectory does not exist!');
    end
    addpath(fullfile(functionsDirectory, '/Downloaded_Functions/'));
                                                % for abf2load.m
    addpath(fullfile(functionsDirectory, '/Brians_Functions/'));
                                                % for identify_channels.m
    addpath(fullfile(functionsDirectory, '/Adams_Functions/'));
                                                % for find_data_files.m
end
                                                
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'combine_sweeps';

% Add required inputs to the Input Parser
% TODO

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SweepNumbers', sweepNumbersDefault);
addParameter(iP, 'DataType', dataTypeDefault, ...
    @(x) any(validatestring(x, possibleDataTypes)));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
sweepNumbers = iP.Results.SweepNumbers;
dataTypeUser = validatestring(iP.Results.DataType, ...
                    [possibleDataTypes, {'auto'}]);
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);
verbose = iP.Results.Verbose;

% TODO: Remove the following line when Verbose is implemented
verbose = true;

% Get file identifier from expLabel
fileIdentifier = strrep(strrep(expLabel, '_IPSC', ''), '_EPSC', '');

% Determine data type, list all .abf, .mat or .txt files from 
%   data subdirectory
[dataType, allDataFiles, nDataFiles, message] = ...
    find_data_files (dataTypeUser, dataDirectory, possibleDataTypes, ...
                        'FileIdentifier', fileIdentifier);

% Check if files are found
if nDataFiles <= 0
    message = {message, 'Failed to combine sweeps!'};
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

% Concatenate sweep data
allData = [];                       % Initialize combined sweep data
if strcmpi(dataMode, 'Katie')       % if there is only one sweep per file
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
        channelTypes = identify_channels(data);
        if isempty(channelTypes{1})
            idxCurrent = 1;
        else
            idxCurrent = strcmp('current', channelTypes);
        end

        % Get current vector (usually pA) for this sweep
        % Note: assume the 2nd dimension is the channel
        current = data(:, idxCurrent);

        % Add to allData
        allData = [allData; current];
    end
elseif strcmpi(dataMode, 'Peter')   % if there are multiple sweeps per file
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
    channelTypes = identify_channels(data);
    idxCurrent = strcmp('current', channelTypes);
  
    % Do for each sweep
    for iSwp = sweepNumbers
        % Note: assume the 2nd dimension is the channel
        %       and the 3rd dimension is the sweep
        current = data(:, idxCurrent, iSwp);

        % Add to allData
        allData = [allData; current];
    end
end

% Show success message
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

%}
