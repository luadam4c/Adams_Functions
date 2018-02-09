function print_or_show_message(toShow, message, varargin)
% Either print a message in standard output or show a message box
% Usage: print_or_show_message(toShow, message, varargin)
% Arguments: TODO: Complete this in the style of function_template.m
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
%   2018-02-02 Created
%   2018-02-07 MD - Input Parser implemented, function working.
%   2018-02-08 AL - Changed specification of toShow from isnumeric
%                       to both isnumeric and islogical
%   2018-02-08 AL - Reverted back to using uiwait() and 'modal'

%% Default values for optional arguments
mTitleDefault = 'Message box';      % default : Title for message box will
                                    %   display as 'Message box'.
iconDefault = 'none';               % default : Does not display an icon with
                                    %   with message box.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = 'print_or_show_message';
iP.KeepUnmatched = false;

% Add required inputs to the Input Parser
addRequired(iP, 'toShow', ...               % whether to show message box
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addRequired(iP, 'message', ...              % the message to show
    @(x) iscellstr(x) || ischar(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MTitle', mTitleDefault, @ischar);
addParameter(iP, 'Icon', iconDefault, @ischar);

% Read from the Input Parser
parse(iP, toShow, message, varargin{:});
mTitle = iP.Results.MTitle;
icon = iP.Results.Icon;

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


%{
OLD CODE:

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

%}
