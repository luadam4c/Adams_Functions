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
%                       'abf'       - AXON binary files
%                       'mat'       - MATLAB data files
%
% Requires:
%       /home/Matlab/Downloaded_Functions/abf2load.m
%       /home/Matlab/Brians_Functions/identify_channels.m
%       /home/Matlab/Adams_Functions/find_data_files.m
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
% 

%% Lists
possibleDataTypes = {'abf', 'mat', 'txt'};     
                    % Precedence: .abf > .mat > .txt

%% Default values for optional arguments
sweepNumbersDefault = 'all';
dataTypeDefault = 'default';    % to detect input data type 
                                %   from possibleDataTypes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsDirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsDirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsDirectory does not exist!');
end
if ~isdeployed
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

% Read from the Input Parser
parse(iP, varargin{:});
sweepNumbers = iP.Results.SweepNumbers;
dataTypeUser = validatestring(iP.Results.DataType, ...
                    [possibleDataTypes, {'default'}]);

% Get file identifier from expLabel
fileIdentifier = strrep(strrep(expLabel, '_IPSC', ''), '_EPSC', '');

% Determine data type, list all .abf, .mat or .txt files from 
%   data subdirectory
[dataType, allDataFiles, nDataFiles, message] = ...
    find_data_files (dataTypeUser, dataDirectory, possibleDataTypes, ...
                        'FileIdentifier', fileIdentifier);
fprintf('%s\n', message);

% Check if files are found
if nDataFiles <= 0
    fprintf('Cannot find data files! Failed to combine sweeps!\n\n');
    return;
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
        fprintf('Too few data files! Failed to combine sweeps!\n\n');
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
        case 'mat'
            % Import data
            data = importdata(dataFileName);
        case 'txt'
            % Import data (the current vector must be one of the columns)
            %   Only one cell per file!
            data = importdata(dataFileName);            
        otherwise
            fprintf(['The data type .%s is not supported yet!\n', ...
                    ' Failed to combine sweeps!\n\n'], dataType);
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
        fprintf('Too many data files! Failed to combine sweeps!\n\n');
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

fprintf('Sweeps successfully combined!!\n\n');

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

%}
