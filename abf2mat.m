function abf2mat(abfFileOrDir, varargin)
%% Converts .abf files to .mat files with time vector (in ms) included
% NOTE: Currently, data must be 2-dimensional (1 sweep per channel) 
%           if saveIndividual is true (default). Otherwise, set
%           'OmitTime' to true and 'SaveIndividual' to false to
%           simply save the data matrix
% Usage: abf2mat(abfFileOrDir, varargin)
%
% Examples: abf2mat('DBA01_31_2016_cage1.abf');
%           abf2mat('Voltage_clamp', ...
%                   'OmitTime', false, ...
%                   'SaveIndividual', true, 'SaveTogether', false, ...
%                   'MaxRecordingLength', 60 * 60 * 1000, ...
%                   'AnimalStr', '_rat');
%
% Arguments:
%       abfFileOrDir- file name of .abf file or 
%                       file directory containing .abf files
%                   if not full path, file must be in current directory
%                   must be a string scalar or a character vector
%       varargin    - 'OmitTime': whether to omit time vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OmitChannels': array of channels to omit
%                   must be a positive integer vector
%                   default == []
%                   - 'SaveIndividual': whether to save individual vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveTogether': whether to save vectors together
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxRecordingLength': maximum recording length in ms
%                   must be a positive scalar
%                   default == Inf
%                   - 'ChannelsPerAnimal': number of channels per animal
%                   must be a positive integer scalar
%                   default == []
%                   - 'AnimalStr': string in file names that separate animals
%                   must be a string scalar or a character vector
%                   default == '_animal'
%                   - 'ChannelStr': string in file names that separate channels
%                   must be a string scalar or a character vector
%                   default == '_animal'
%                   - 'PieceStr': string in file names that separate pieces
%                   must be a string scalar or a character vector
%                   default == '_animal'
%                   - 'OutFolder': directory to output csv file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == 'matfiles'
% Requires:
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%       /home/Matlab/Downloaded_Functions/abf2load.m
%       cd/check_dir.m


% File History:
% 2016-08-02 Created
% 2016-09-19 Modified so that it reads files in 
%               current directory without full path
% 2016-09-19 Modified so that it outputs files 
%               containing each channel individually as well
% 2017-08-01 Changed tabs to spaces
% 2017-08-01 Added the option of not containing time vectors
% 2017-08-01 Added input parser scheme
% 2017-08-01 Made saving individual vectors and option
% 2017-08-01 Now accepts a directory as the first argument
% 2018-04-26 Add 'MaxRecordingLength' as a parameter
% 2018-04-30 Changed num2str(iPiece) -> num2str(iPiece, '%02.f')
%               so that files will be in alphabetical order
% 2018-05-01 Changed num2str(iChannel) -> num2str(iChannel, '%02.f')
%               so that files will be in alphabetical order
% 2018-05-02 Added outFolder and made the default 'matfiles'
% 2018-05-17 Added ChannelsPerAnimal and animalStr
% 2018-05-17 Added OmitChannels
% TODO: Apply identify_channels.m and save each channel
% TODO: Update this to use parse_abf.m

%% Hard-coded parameters

%% Default values for optional arguments
omitTimeDefault = false;        % whether to omit time vector by default
omitChannelsDefault = [];       % list of channels to omit by default
saveIndividualDefault = true;   % whether to save individual vectors by default
saveTogetherDefault = false;    % whether to save vectors together by default
maxRecordingLengthDefault = Inf;
channelsPerAnimalDefault = [];  % default number of channels per animal
animalStrDefault = '_animal';   % string in file names that separate animals
channelStrDefault = '_channel'; % string in file names that separate channels
pieceStrDefault = '_piece';     % string in file names that separate pieces
outFolderDefault = 'matfiles';  % default directory to output mat file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions across servers
if ~isdeployed
    if exist('/home/Matlab/', 'dir') == 7
        functionsDirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsDirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsDirectory does not exist!');
    end
    addpath(fullfile(functionsDirectory, 'Downloaded_Functions')); 
                                            % for abf2load.m
end

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help abf2mat'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'abfFileOrDir', ... % file name or file directory
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OmitTime', omitTimeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OmitChannels', omitChannelsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'vector'}));
addParameter(iP, 'SaveIndividual', saveIndividualDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveTogether', saveTogetherDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxRecordingLength', maxRecordingLengthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'ChannelsPerAnimal', channelsPerAnimalDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'AnimalStr', animalStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ChannelStr', channelStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PieceStr', pieceStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, abfFileOrDir, varargin{:});
omitTime = iP.Results.OmitTime;
omitChannels = iP.Results.OmitChannels;
saveIndividual = iP.Results.SaveIndividual;
saveTogether = iP.Results.SaveTogether;
maxRecordingLength = iP.Results.MaxRecordingLength;
channelsPerAnimal = iP.Results.ChannelsPerAnimal;
animalStr = iP.Results.AnimalStr;
channelStr = iP.Results.ChannelStr;
pieceStr = iP.Results.PieceStr;
outFolder = iP.Results.OutFolder;

% Parse first argument
if exist(abfFileOrDir, 'file') == 2                         % it's a file
    % Store the directory containing the file
    [abfDir, abfFileBase, abfFileExt] = fileparts(abfFileOrDir);

    % Get the relative path for the file name
    abfFileName = [abfFileBase, abfFileExt];

    % Set flag
    multipleFiles = false;
elseif exist(abfFileOrDir, 'dir') == 7                      % it's a directory
    % The argument is already the full directory
    abfDir = abfFileOrDir;

    % Set flag
    multipleFiles = true;
elseif exist(fullfile(pwd, abfFileOrDir), 'file') == 2      % it's a file
    % The first argument is just the file name in current directory
    abfDir = pwd;

    % Get the relative path for the file name
    abfFileName = abfFileOrDir;

    % Set flag
    multipleFiles = false;
elseif exist(fullfile(pwd, abfFileOrDir), 'dir') == 7       % it's a directory
    % The first argument is a subdirectory in current directory
    %   Get full path to directory
    abfDir = fullfile(pwd, abfFileOrDir);

    % Set flag
    multipleFiles = true;
else
    message = sprintf('The .abf file or directory %s does not exist!', ...
                        abfFileOrDir);
    mTitle = 'File or Directory Not Found';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                          'MessageMode', 'show', 'Verbose', true, ...
                          'CreateMode', 'replace');
    return;
end

% Check if output directory exists
check_dir(outFolder);

% Find all .abf files to convert
if multipleFiles
    % List all the .abf files in the directory
    allAbfFiles = dir(fullfile(abfDir, '*.abf'));

    % Compute the number of .abf files in the directory
    nAbfFiles = length(allAbfFiles);

    % Loop through all files
%    for iFile = 1:nAbfFiles
    parfor iFile = 1:nAbfFiles
        abfFileName = allAbfFiles(iFile).name;
        convert_abf2mat(abfDir, abfFileName, outFolder, ...
                        omitTime, omitChannels, ...
                        saveIndividual, saveTogether, ...
                        maxRecordingLength, channelsPerAnimal, ...
                        animalStr, channelStr, pieceStr);
    end
else
    convert_abf2mat(abfDir, abfFileName, outFolder, ...
                    omitTime, omitChannels, ...
                    saveIndividual, saveTogether, ...
                    maxRecordingLength, channelsPerAnimal, ...
                    animalStr, channelStr, pieceStr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_abf2mat(abfDir, abfFileName, outFolder, omitTime, omitChannels, saveIndividual, saveTogether, maxRecordingLength, channelsPerAnimal, animalStr, channelStr, pieceStr)
%% Convert an .abf file to matfile(s)
% TODO: work with 3-D data

% Load abf file
%   d       - the data with 3 voltage channels; 
%   siUs    - the sampling interval in us
abfFullFileName = fullfile(abfDir, abfFileName);
[d, siUs, ~] = abf2load(abfFullFileName);

nDps = size(d, 1);          % Total number of data points in each channel
nChannelsAll = size(d, 2);  % Total number of channels in each abf file
siMs = siUs/1000;           % sampling interval in ms

% Compute the total recording length in ms
totalRecordingLength = siMs * nDps;

% Decide if there is more than one animal
if ~isempty(channelsPerAnimal) && nChannelsAll > channelsPerAnimal
    moreThanOneAnimal = true;

    % Compute the number of animals
    nAnimals = floor(nChannelsAll / channelsPerAnimal);
else
    moreThanOneAnimal = false;
    nAnimals = 1;
    channelsPerAnimal = nChannelsAll;
end

% Decide whether to break up the trace into multiple files:
%   If greater than max recording length, break up the trace into pieces
if totalRecordingLength > maxRecordingLength
    breakUpTrace = true;

    % Compute the number of pieces after breaking up the trace
    nPieces = ceil(totalRecordingLength / maxRecordingLength);

    % Compute the number of indices per piece
    indPerPiece = round(maxRecordingLength / siMs);
else
    breakUpTrace = false;
    nPieces = 1;
    indPerPiece = nDps;
end

% For each piece
for iPiece = 1:nPieces
    % Get the index start and end points for this piece
    indThisStart = (iPiece - 1) * indPerPiece + 1;
    if iPiece < nPieces
        indThisEnd = (iPiece - 1) * indPerPiece + indPerPiece;
    elseif iPiece == nPieces
        indThisEnd = nDps;
    else
        error('error with code!');
    end

    % Construct the full piece string if any
    if breakUpTrace
        fullPieceStr = [pieceStr, num2str(iPiece, '%02.f')];
    else
        fullPieceStr = '';
    end        

    % Get the time vector in ms
    t = siMs*(indThisStart:indThisEnd)';

    % For each animal
    for iAnimal = 1:nAnimals
        % Construct the full animal string if any
        if moreThanOneAnimal
            fullAnimalStr = [animalStr, num2str(iAnimal)];
        else
            fullAnimalStr = '';
        end

        % Find the channels that the user wants to save
        channelsToSave = setdiff(1:channelsPerAnimal, omitChannels);

        % Find the original channel numbers in the abf file
        if moreThanOneAnimal
            channelNos = (iAnimal - 1) * channelsPerAnimal + channelsToSave;
        else
            channelNos = channelsToSave;
        end

        % Get the number of distinct channel numbers
        nChannelsToSave = length(channelNos);

        % Save individual vectors into its own mat file (default)
        if saveIndividual
            for iChannel = 1:nChannelsToSave
                % Get the original channel number in the abf file
                channelNo = channelNos(iChannel);

                % Get the channel number in the animal
                channelToSave = channelsToSave(iChannel);

                % Save as 'vec'
                vec = d(indThisStart:indThisEnd, channelNo);

                % Construct the full channel string
                fullChannelStr = [channelStr, num2str(channelToSave, '%02.f')];

                % Create the matfile name
                subMatFileName = strrep(abfFileName, '.abf', ...
                                        [fullAnimalStr, fullChannelStr, ...
                                        fullPieceStr, '.mat']);                

                % Save vectors to matfile
                if ~omitTime
                    save(fullfile(outFolder, subMatFileName), ...
                         't', 'vec', '-v7.3');
                else
                    save(fullfile(outFolder, subMatFileName), 'vec');
                end
            end
        end

        % Create a data vector for each channel
        for iChannel = 1:nChannelsToSave
            % Get the original channel number in the abf file
            channelNo = channelNos(iChannel);

            % Get the channel number in the animal
            channelToSave = channelsToSave(iChannel);

            % Name by the channel number for this animal
            %   (in uV for EEG)
            eval(sprintf('vec%d = d(%d:%d, %d);', ...
                        channelToSave, indThisStart, indThisEnd, channelNo));
        end

        % Save all vectors into one v7.3 mat file
        if saveTogether
            % Create matfile name
            matFileName = strrep(abfFileName, '.abf', ...
                                [fullAnimalStr, fullPieceStr, '.mat']);

            % Save vectors to matfile
            command_string_part = [];
            for iChannel = 1:nChannelsToSave
                % Get the channel number in the animal
                channelToSave = channelsToSave(iChannel);

                % Construct the part of the command string that saves 
                %   this channel
                command_string_part = [command_string_part, ...
                                        ', ''vec', num2str(channelToSave), ''''];
            end
            if ~omitTime
                eval(sprintf(['save(fullfile(outFolder, matFileName)', ...
                                ', ''t''%s, ''-v7.3'');'], ...
                                command_string_part));
            else
                eval(sprintf('save(fullfile(outFolder, matFileName), ''d'');'));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

v1 = d(:, 1);                % voltage vector, channel #1
v2 = d(:, 2);                % voltage vector, channel #2
v3 = d(:, 3);                % voltage vector, channel #3

save(matFullFileName, 't', 'v1', 'v2', 'v3', '-v7.3');

% nChannels = size(d, 2);     % Total number of channels in each abf file

% Create a matfile name with the piece number
if breakUpTrace
else
    % Create a matfile name
    matFileName = strrep(abfFileName, '.abf', '.mat');
end

if breakUpTrace
    % Create a matfile name with the channel & piece number
else
    % Create a matfile name with the channel number
    subMatFileName = strrep(abfFileName, '.abf', ...
                    [fullChannelStr, '.mat']);
end

% If more than one animal, get the actual channel number
if moreThanOneAnimal
    iChannel = mod(channelNo, channelsPerAnimal);
else
    iChannel = channelNo;
end

% Compute the number of channels that are for animals
nChannels = nAnimals * channelsPerAnimal;
nChannels = nChannelsAll;

% The argument is already the full path
abfFullFileName = abfFileOrDir;
    %   Get full path to file and directory
abfFullFileName = fullfile(pwd, abfFileOrDir);

% Set dependent argument defaults
if isempty(outFolder)
    % Default output directory is matfiles under the data directory
    outFolder = fullfile(abfDir, 'matfiles');
end

% Make sure the output directory exists
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%