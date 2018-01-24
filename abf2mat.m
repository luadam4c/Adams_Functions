function abf2mat(abfFileOrDir, varargin)
%% Converts .abf files to .mat files with time vector included
% NOTE: Currently, data must be 2-dimensional (1 sweep per channel) 
%           is saveIndividual is true (default). Otherwise, set
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
% TODO: Apply identify_channels.m and save each channel

%% Default values for optional arguments
omitTimeDefault = false;        % whether to omit time vector by default
saveIndividualDefault = true;   % whether to save individual vectors by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help abf2mat'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'abf2mat';

% Add required inputs to the Input Parser
addRequired(iP, 'abfFileOrDir', ... % file name or file directory
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OmitTime', omitTimeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveIndividual', saveIndividualDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, abfFileOrDir, varargin{:});
omitTime = iP.Results.OmitTime;
saveIndividual = iP.Results.SaveIndividual;

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

% Find all .abf files to convert
if multipleFiles
    allAbfFiles = dir(fullfile(abfDir, '*.abf'));
    nAbfFiles = length(allAbfFiles);   % # of .abf files in directory
    for iFile = 1:nAbfFiles
        convert_abf2mat(fullfile(abfDir, allAbfFiles(iFile).name), ...
                        omitTime, saveIndividual);
    end
else
    convert_abf2mat(abfFullFileName, omitTime, saveIndividual);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_abf2mat(abfFullFileName, omitTime, saveIndividual)
%% TODO: work with 3-D data

% Load abf file
%   d       - the data with 3 voltage channels; 
%   siUs    - the sampling interval in us
[d, siUs, ~] = abf2load(abfFullFileName);    

nDps = size(d, 1);          % Total number of data points in each channel
nChannels = size(d, 2);     % Total number of channels in each abf file
siMs = siUs/1000;           % sampling interval in ms
t = siMs*(1:nDps)';         % time vector in ms

% Create a data vector 
for k = 1:nChannels
    eval(sprintf('vec%d = d(:, %d);', k, k));   % voltage vector 
                                                %   (in uV for EEG), channel #k
end

% Save all vectors into one v7.3 mat file
matFullFileName = strrep(abfFullFileName, '.abf', '.mat');
command_string_part = [];
for k = 1:nChannels
    command_string_part = [command_string_part, ', ''vec', num2str(k), ''''];
end
if ~omitTime
    eval(sprintf('save(matFullFileName, ''t''%s, ''-v7.3'');', ...
                    command_string_part));
else
    eval(sprintf('save(matFullFileName, ''d'');'));
end

% Save individual vectors into its own mat file (old version)
if saveIndividual
    for k = 1:nChannels
        vec = d(:, k);
        subMatFullFileName = strrep(abfFullFileName, '.abf', ...
                                ['_channel', num2str(k), '.mat']);
        if ~omitTime
            save(subMatFullFileName, 't', 'vec', '-v7.3');
        else
            save(subMatFullFileName, 'vec');
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
