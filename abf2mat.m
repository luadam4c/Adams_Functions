function abf2mat(abfFileOrDir, varargin)
%% Converts .abf files to .mat files with time vector (in ms) included
% NOTE: Currently, data must be 2-dimensional (1 sweep per channel) 
%           if saveIndividual is true (default). Otherwise, set
%           'OmitTime' to true and 'SaveIndividual' to false to
%           simply save the data matrix
% Usage: abf2mat(abfFileOrDir, varargin)
%
% Examples: abf2mat('DBA01_31_2016_cage1.abf')
%           abf2mat('Voltage_clamp', 'OmitTime', true, 'SaveIndividual', false)
%
% Arguments:
%       abfFileOrDir- file name of .abf file or 
%                       file directory containing .abf files
%                   if not full path, file must be in current directory
%                   must be a string scalar or a character vector
%       varargin    - 'OmitTime': whether to omit time vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveIndividual': whether to save individual vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveTogether': whether to save vectors together
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxRecordingLength': maximum recording length in ms
%                   must be a positive number
%                   default == Inf
%                   - 'OutFolder': directory to output csv file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == 'fullfile(abfDir, matfiles)'

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
% TODO: Apply identify_channels.m and save each channel

%% Hard-coded parameters
% TODO: Make these default values for optional arguments
channelStr = '_channel';        % string in file names that separate channels
pieceStr = '_piece';            % string in file names that separate pieces

%% Default values for optional arguments
omitTimeDefault = false;        % whether to omit time vector by default
saveIndividualDefault = true;   % whether to save individual vectors by default
saveTogetherDefault = false;    % whether to save vectors together by default
maxRecordingLengthDefault = Inf;
outFolderDefault = '';          % default directory to output csv file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
addParameter(iP, 'SaveIndividual', saveIndividualDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveTogether', saveTogetherDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxRecordingLength', maxRecordingLengthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, abfFileOrDir, varargin{:});
omitTime = iP.Results.OmitTime;
saveIndividual = iP.Results.SaveIndividual;
saveTogether = iP.Results.SaveTogether;
maxRecordingLength = iP.Results.MaxRecordingLength;
outFolder = iP.Results.OutFolder;

% Parse first argument
if exist(abfFileOrDir, 'file') == 2                         % it's a file
    % The argument is already the full path
    abfFullFileName = abfFileOrDir;

    % Store directory containing the file
    [abfDir, ~, ~] = fileparts(abfFileOrDir);

    % Set flag
    multipleFiles = false;
elseif exist(abfFileOrDir, 'dir') == 7                      % it's a directory
    % The argument is already the full directory
    abfDir = abfFileOrDir;

    % Set flag
    multipleFiles = true;
elseif exist(fullfile(pwd, abfFileOrDir), 'file') == 2      % it's a file
    % The first argument is just the file name in current directory
    %   Get full path to file and directory
    abfFullFileName = fullfile(pwd, abfFileOrDir);
    abfDir = pwd;

    % Set flag
    multipleFiles = false;
elseif exist(fullfile(pwd, abfFileOrDir), 'dir') == 7       % it's a directory
    % The first argument is a subdirectory in current directory
    %   Get full path to directory
    abfDir = fullfile(pwd, abfFileOrDir);

    % Set flag
    multipleFiles = true;
else
    error('File or directory undefined!');
end

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

% Find all .abf files to convert
if multipleFiles
    allAbfFiles = dir(fullfile(abfDir, '*.abf'));
    nAbfFiles = length(allAbfFiles);   % # of .abf files in directory
    for iFile = 1:nAbfFiles
        abfFileName = allAbfFiles(iFile).name;
        convert_abf2mat(abfDir, abfFileName, outFolder, ...
                        omitTime, saveIndividual, saveTogether, ...
                        maxRecordingLength);
    end
else
    convert_abf2mat(abfDir, abfFileName, outFolder, ...
                    omitTime, saveIndividual, saveTogether, ...
                    maxRecordingLength);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_abf2mat(abfDir, abfFileName, outFolder, omitTime, saveIndividual, saveTogether, maxRecordingLength)
%% TODO: work with 3-D data

% Load abf file
%   d       - the data with 3 voltage channels; 
%   siUs    - the sampling interval in us
abfFullFileName = fullfile(abfDir, abfFileName);
[d, siUs, ~] = abf2load(abfFullFileName);

nDps = size(d, 1);          % Total number of data points in each channel
nChannels = size(d, 2);     % Total number of channels in each abf file
siMs = siUs/1000;           % sampling interval in ms

% Compute the total recording length in ms
totalRecordingLength = siMs * nDps;

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

    % Get the time vector in ms
    t = siMs*(indThisStart:indThisEnd)';

    % Save individual vectors into its own mat file (default)
    if saveIndividual
        for iChannel = 1:nChannels
            % Save as 'vec'
            vec = d(indThisStart:indThisEnd, iChannel);

            if breakUpTrace
                % Create a matfile name with the channel & piece number
                subMatFileName = strrep(abfFileName, '.abf', ...
                            [channelStr, num2str(iChannel, '%02.f'), ...
                                pieceStr, num2str(iPiece, '%02.f'), '.mat']);                
            else
                % Create a matfile name with the channel number
                subMatFileName = strrep(abfFileName, '.abf', ...
                            [channelStr, num2str(iChannel, '%02.f'), '.mat']);
            end
            if ~omitTime
                save(fullfile(outFolder, subMatFileName), 't', 'vec', '-v7.3');
            else
                save(fullfile(outFolder, subMatFileName), 'vec');
            end
        end
    end

    % Create a data vector for each channel
    for iChannel = 1:nChannels
        % Name by the channel number
        %   (in uV for EEG), channel #iChannel
        eval(sprintf('vec%d = d(%d:%d, %d);', ...
                    iChannel, indThisStart, indThisEnd, iChannel));
    end

    % Save all vectors into one v7.3 mat file
    if saveTogether
        if breakUpTrace
            % Create a matfile name with the piece number
            matFileName = strrep(abfFileName, '.abf', ...
                                [pieceStr, num2str(iPiece, '%02.f'), '.mat']);
        else
            % Create a matfile name
            matFileName = strrep(abfFileName, '.abf', '.mat');
        end
        command_string_part = [];
        for k = 1:nChannels
            command_string_part = [command_string_part, ...
                                    ', ''vec', num2str(k), ''''];
        end
        if ~omitTime
            eval(sprintf('save(fullfile(outFolder, matFileName), ', ...
                            '''t''%s, ''-v7.3'');', ...
                            command_string_part));
        else
            eval(sprintf('save(fullfile(outFolder, matFileName), ''d'');'));
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

%}
