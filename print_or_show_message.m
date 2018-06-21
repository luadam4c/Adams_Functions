function print_or_show_message(message, varargin)
%% Either print a message in standard output or show a message box
% Usage: print_or_show_message(message, varargin)
% Explanation: 
%       Either pause program and show message box, only show message box, 
%           or printed in stardard output
% Example:
%       print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
%                             'MessageMode', 'show', 'Verbose', true, ...
%                             'CreateMode', 'non-modal')
% Arguments:
%       message     - message displayed in message box
%                   must be a string scalar or a character vector
%       varargin    - 'MTitle': Title of message box
%                   must be a character vector
%                   default == 'Message box'
%                   - 'Icon': displayed icon on message box
%                   must be a character vector
%                   default == 'none'
%                   - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%                   - 'CreateMode' - window mode for the msgbox() function 
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'non-modal' - does not disable main content and 
%                                       does not overwrite other dialog boxes
%                       'modal'     - disables main content and overwrites
%                                       dialog boxes of the same title
%                       'replace'   - overwrites dialog boxes of the same title
%                                       with non-modal dialog boxes
%                   default == 'modal'
%                   - 'Verbose' - whether to print to standard output
%                                   regardless of message mode
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   
% 
% Requires:
%       /home/Matlab/Miras_Functions/print_cellstr.m
%
% Used by:
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/minEASE/combine_eventInfo.m
%       /home/Matlab/minEASE/compute_plot_average_PSC_traces.m
%       /home/Matlab/minEASE/detect_gapfree_events.m
%       /home/Matlab/Kojis_Functions/compute_average_PSC_trace.m
%       /home/Matlab/Adams_Functions/abf2mat.m
%       /home/Matlab/Adams_Functions/atf2sheet.m
%       /home/Matlab/Adams_Functions/combine_sweeps.m
%       /home/Matlab/Adams_Functions/create_input_file.m
%       /home/Matlab/EEG_gui/combine_EEG_gui_outputs.m
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m
%
% File History:
% 2018-02-02 Created
% 2018-02-07 MD - Input Parser implemented, function working.
% 2018-02-08 AL - Changed specification of toShow from isnumeric
%                   to both isnumeric and islogical
% 2018-02-08 AL - Reverted back to using uiwait() and 'modal'
% 2018-02-26 MD - Modified to have messageMode and removed toShow, messagebox
%                   always shows now with the fprintf being optional.
% 2018-03-01 MD - Added verbose parameter
% 2018-05-15 AL - Added example usage
% 2018-05-15 AL - Added CreateMode

%% Hard-coded parameters
validMessageModes = {'wait', 'show', 'none'};
validCreateModes = {'non-modal', 'modal', 'replace'};

%% Default values for optional arguments
mTitleDefault = 'Message box';      % default : Title for message box will
                                    %   display as 'Message box'.
iconDefault = 'none';               % default : Does not display an icon with
                                    %   with message box.
messageModeDefault = 'wait';        % default : Pauses program and displays
                                    %   message box.
createModeDefault = 'modal';        % default : Disables main content and 
                                    %   overwrites dialog boxes of the same title
verboseDefault = false;             % default: Program does not print message
                                    %   even if message box is shown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions across servers
if ~isdeployed
    if exist(fullfile(pwd, 'Miras_Functions'), 'dir') == 7
        functionsdirectory = pwd;
    elseif exist('/home/Matlab/', 'dir') == 7
        functionsDirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsDirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsDirectory does not exist!');
    end
    addpath(fullfile(functionsDirectory, 'Miras_Functions')); 
                                            % for print_cellstr.m
end

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'message', ...              % the message to show
    @(x) iscellstr(x) || ischar(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MTitle', mTitleDefault, @ischar);
addParameter(iP, 'Icon', iconDefault, @ischar);
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'CreateMode', createModeDefault, ...
    @(x) any(validatestring(x, validCreateModes)));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, message, varargin{:});
mTitle = iP.Results.MTitle;
icon = iP.Results.Icon;
verbose = iP.Results.Verbose;

% Match possibly ambiguous strings to valid strings
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);
createMode = validatestring(iP.Results.CreateMode, validCreateModes);

% Print to standard output if verbose is true or messageMode == 'none'
if verbose || strcmp(messageMode, 'none')
    messageStr = print_cellstr(message, 'Delimiter', '\n', ...
                                        'OmitQuotes', true, ...
                                        'OmitBraces', true, ...
                                        'OmitNewline', false, ...
                                        'ToPrint', true);
end

% Display message box and/or stop program
switch messageMode
case 'wait'
    %   Program stops and displays message box
    uiwait(msgbox(message, mTitle, icon, createMode));                  
case 'show'
    %   Program does not stop but still displays message box
    msgbox(message, mTitle, icon, createMode);
case 'none'
    %   Program does not stop or display message box
otherwise
    error(['This is not a recognized message mode. ', ...
            'Valid message modes are ''wait'', ''show'', or ''none''.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) any(validatestring(x, {'true', 'false'})));

addRequired(iP, 'toShow', ...               % whether to show message box
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Show a message box or print to standard output
if toShow               % if user wants to show a message box
    % Display a message box (replaceable by another box with the same mTitle)
    %   and wait for the user to close it
    uiwait(msgbox(message, mTitle, icon, 'modal'));
else                    % if user does not want to show a message box
    % Print the message to standard output
    messageStr = print_cellstr(message, 'Delimiter', '\n', 'OmitNewline', true);
    fprintf('%s\n', messageStr);
end

    if iscell(message)
        messageStr = strjoin(message, '\n');
    else
        messageStr = message;
    end
    
if toShow
    uiwait(msgbox(message, mTitle, icon, 'modal'));
else
    messageStr = print_cellstr(message, 'Delimiter', '\n', 'OmitNewline', true);
    fprintf('%s\n', messageStr);
end

%Read from Input Parser:
parse(iP, toShow, message, mTitle, icon, varargin{:});
mtitle = iP.Results.mTitle;
Icon = iP.Results.icon;

function print_or_show_message(toShow, message, mTitle, icon, varargin)

% User Input Specifications:
if toShow == 1

    msgbox(message, mTitle, icon);

if verbose && (strcmp(messageMode, 'wait') ||  strcmp(messageMode, 'show'))

case 'none'
    %   Program does not stop or display message box and prints message in
    %   standard output
    messageStr = print_cellstr(message, 'Delimiter', '\n', 'OmitNewline', true);
    fprintf('%s\n', messageStr); 

else
% warning('off', verbose) isn't this how verbose is used? I also looked in
% minEASE to see how it was used but I'm still confused. I was able to add it
% as a pair parameter though!

uiwait(msgbox(message, mTitle, icon, 'modal'));                  
msgbox(message, mTitle, icon, 'modal');

messageStr = print_cellstr(message, 'Delimiter', '\n', ...
                                    'OmitQuotes', true, ...
                                    'OmitBraces', true, ...
                                    'OmitNewline', true, ...
                                    'ToPrint', false);
fprintf('%s\n', messageStr);
%}
