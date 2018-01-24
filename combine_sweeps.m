function allData = combine_sweeps(dataDirectory, expLabel, dataMode, varargin)
%% Combines sweeps that begin with expLabel in dataDirectory under dataMode
% Usage: allData = combine_sweeps(dataDirectory, expLabel, dataMode, varargin)
% Explanation:
%       TODO
% Outputs:
%       TODO
% Arguments:    
%       TODO
%
% Requires:
%       /home/Matlab/Downloaded_Functions/abf2load.m
%       /home/Matlab/Brians_Functions/identify_channels.m
%
% Used by:    
%       /home/Matlab/minEASE/minEASE.m
%
% File History:
% 2017-07-25 Created by AL
% 2017-10-15 Added success message
% 

%% Default values for optional arguments
sweepNumbersDefault = 'all';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
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
                                                
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'combine_sweeps';

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SweepNumbers', sweepNumbersDefault);

% Read from the Input Parser
parse(iP, varargin{:});
sweepNumbers = iP.Results.SweepNumbers;

% Get file identifier from expLabel
fileIdentifier = strrep(strrep(expLabel, '_IPSC', ''), '_EPSC', '');

% Find all the .abf files in the dataDirectory with this fileIdentifier
abfFiles = dir(fullfile(dataDirectory, [fileIdentifier, '*.abf']));
nAbfFiles = length(abfFiles);

% Check if files are found
if nAbfFiles <= 0 
    fprintf('Cannot find abf files! Failed to combine sweeps!\n\n');
    return;
end

% Concatenate sweep data
allData = [];                       % Initialize combined sweep data
if strcmpi(dataMode, 'Katie')       % if there is only one sweep per file
    % Set default sweep numbers for one sweep per file
    if ~isnumeric(sweepNumbers) && strcmpi(sweepNumbers, 'all')
        sweepNumbers = 1:nAbfFiles;
    end

    % Check if number of abf files is correct
    if nAbfFiles < max(sweepNumbers)
        fprintf('Too few abf files! Failed to combine sweeps!\n\n');
        return;
    end

    % Do for each sweep/file
    for iFile = sweepNumbers
        % Load abf data for this sweep
        abfFileName = fullfile(dataDirectory, abfFiles(iFile).name);
        [data, ~] = abf2load(abfFileName);

        % Identify the channels in data
        channelTypes = identify_channels(data);
        idxCurrent = strcmp('current', channelTypes);

        % Get current vector (usually pA) for this sweep
        % Note: assume the 2nd dimension is the channel
        current = data(:, idxCurrent);

        % Add to allData
        allData = [allData; current];
    end
elseif strcmpi(dataMode, 'Peter')   % if there are multiple sweeps per file
    % Check if multiple abf files exist
    if nAbfFiles > 1
        fprintf('Too many abf files! Failed to combine sweeps!\n\n');
        return;
    end

    % Load abf data for all sweeps
    abfFileName = fullfile(dataDirectory, abfFiles(1).name);
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

fprintf('Sweeps successfully combined!!\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

