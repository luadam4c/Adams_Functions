function print_and_show_message(message, varargin)
%% Print to standard output and show message box at the same time
% Usage: print_and_show_message(message, varargin)
% Explanation: 
%       Either pause program and show message box, only show message box, 
%               or neither
%       Message will be printed in stardard output in all cases.
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
%       /home/Matlab/Adams_Functions/combine_sweeps.m
%
% File History:
%   2018-02-14 Modified from print_or_show_message.m
%   2018-02-14 MD - Added optional cases for messageMode
%   2018-02-14 MD - Improved documentation
%   2018-02-14 AL - Changed description and usage documentation

%% Hard-coded parameters
validMessageModes = {'wait', 'show', 'none'};

%% Default values for optional arguments
mTitleDefault = 'Message box';      % default : Title for message box will
                                    %   display as 'Message box'.
iconDefault = 'none';               % default : Does not display an icon with
                                    %   with message box
messageModeDefault = 'wait';        % default : Pauses program and displays
                                    %   message box

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
    addpath(fullfile(functionsDirectory, 'Miras_Functions')); 
                                            % for print_cellstr.m
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = false;

% Add required inputs to the Input Parser
addRequired(iP, 'message', ...              % the message to show
    @(x) iscellstr(x) || ischar(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MTitle', mTitleDefault, @ischar);
addParameter(iP, 'Icon', iconDefault, @ischar);
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));

% Read from the Input Parser
parse(iP, message, varargin{:});
mTitle = iP.Results.MTitle;
icon = iP.Results.Icon;

% Match possibly ambiguous strings to valid strings
messageMode = validatestring(iP.Results.MessageMode, validMessageModes); 

% Message prints in standard output
messageStr = print_cellstr(message, 'Delimiter', '\n', 'OmitNewline', true);
fprintf('%s\n', messageStr); 
    
%   Display message box and/or stop program
switch messageMode
case 'wait'
    %   Program stops and displays message box               
    uiwait(msgbox(message, mTitle, icon, 'modal'));                  
case 'show'
    %   Program does not stop but still displays message box
    msgbox(message, mTitle, icon, 'modal');
case 'none'
    %   Program does not stop or display message box
otherwise
    error(['This is not a recognized message mode. ', ...
            'Valid message modes are ''wait'', ''show'', or ''none''.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

addRequired(iP, 'messageMode', ...              
    @(x) validateattributes(x, {'char'}, {'nonempty'}));

if messageMode               % if user wants to show a message box
    % Display a message box (replaceable by another box with the same mTitle)
    %   and wait for the user to close it
    uiwait(msgbox(message, mTitle, icon, 'modal'));
else                    % if user does not want to show a message box
    % Print the message to standard output
    messageStr = print_cellstr(message, 'Delimiter', '\n', 'OmitNewline', true);
    fprintf('%s\n', messageStr);
end


addRequired(iP, 'messageMode', ...              
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

    if iscell(message)
        messageStr = strjoin(message, '\n');
    else
        messageStr = message;
    end
    
if messageMode
    uiwait(msgbox(message, mTitle, icon, 'modal'));
else
    messageStr = print_cellstr(message, 'Delimiter', '\n', 'OmitNewline', true);
    fprintf('%s\n', messageStr);
end

%Read from Input Parser:
parse(iP, messageMode, message, mTitle, icon, varargin{:});
mtitle = iP.Results.mTitle;
Icon = iP.Results.icon;

function print_or_show_message(messageMode, message, mTitle, icon, varargin)

% User Input Specifications:
if messageMode == 1

    msgbox(message, mTitle, icon);

%}
