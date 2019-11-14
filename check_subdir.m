function check_subdir (parentDirectory, subDirectories, varargin)
%% Checks if needed subdirectory(ies) exist in parentDirectory
% Usage: check_subdir (parentDirectory, subDirectories, varargin)
% Arguments:
%       parentDirectory - parent directory of subdirectories
%                       must be a scalar text
%       subDirectories - directory(ies) to check
%                       must be a cell array of character arrays
%                           or a scalar text
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%
% Requires:
%       cd/check_dir.m
%
% Used by:    
%       cd/create_subdir_copy_files.m
%       cd/create_input_file.m
%       cd/find_passive_params.m
%       cd/find_istart.m
%       cd/find_IPSC_peak.m
%       cd/find_LTS.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m
%       cd/parse_multiunit.m
%       cd/save_all_zooms.m
%       /media/adamX/m3ha/data_dclamp/take4/find_special_cases.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting22.m
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/EEG_gui/EEG_gui.m
%       /home/Matlab/function_template.m

% File History:
% 2016-11-02 Created
% 2018-06-19 Changed tabs to spaces
% 2018-09-18 Added input parser and verbose
% 2018-10-03 Added 'MessageMode'

%% Hard-coded parameters
validMessageModes = {'wait', 'show', 'none'};

%% Default values for optional arguments
verboseDefault = true;
messageModeDefault = 'none';        % do not display message box by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'parentDirectory', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addRequired(iP, 'subDirectories', ...
    @(x) iscellstr(subDirectories) || ischar(subDirectories) || ...
        isstring(subDirectories));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));

% Read from the Input Parser
parse(iP, parentDirectory, subDirectories, varargin{:});
verbose = iP.Results.Verbose;
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);

%% Construct the full paths for checking
if iscell(subDirectories)
    directories = ...
        cellfun(@(x) fullfile(parentDirectory, x), subDirectories, ...
                'UniformOutput', false);
else
    % Construct the full paths for checking
    directories = fullfile(parentDirectory, subDirectories);
end

%% Check directory(ies)
check_dir(directories, 'Verbose', verbose, 'MessageMode', messageMode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscell(subDirectories)
    for k = 1:numel(subDirectories)
        if exist(fullfile(parentDirectory, subDirectories{k}), 'dir') ~= 7
            mkdir(fullfile(parentDirectory, subDirectories{k}));
            if verbose
                fprintf('New subdirectory is made: %s\n\n', ...
                    fullfile(parentDirectory, subDirectories{k}));
            end
        end
    end
else
    if exist(fullfile(parentDirectory, subDirectories), 'dir') ~= 7
        mkdir(fullfile(parentDirectory, subDirectories));
        if verbose
            fprintf('New subdirectory is made: %s\n\n', ...
                    fullfile(parentDirectory, subDirectories));
        end
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
