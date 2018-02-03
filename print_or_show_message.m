function print_or_show_message(toShow, message, mTitle, icon, varargin)
% Either print a message in standard output or show a message box
% Usage: print_or_show_message(toShow, message, mTitle, icon, varargin)
% 
% Requires:
% TODO:      /home/Matlab/Miras_Functions/print_cellstr.m
%
% Used by:
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/minEASE/combine_eventInfo.m, 
%       /home/Matlab/minEASE/compute_plot_average_PSC_traces.m
%       /home/Matlab/Kojis_Functions/compute_average_PSC_trace.m
%       /home/Matlab/Adams_Functions/combine_sweeps.m
%
% File History:
%   2018-02-02 Created
%   TODO: Set up input Parser
%   TODO: Make mTitle and icon optional arguments
%       with parameter names 'MessageTitle' and 'Icon' and
%       with default values 'Message box' and 'none', respectively
%

%% Default values for optional arguments
mTitleDefault = 'Message box';      % default TODO: Description of MessageTitle
iconDefault = 'none';               % default TODO: Description of Icon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Do the job
if toShow
    uiwait(msgbox(message, mTitle, icon, 'modal'));
else
% TODO: print_cellstr(message, 'Delimiter', '\n', 'OmitNewline', false);
    if iscell(message)
        messageStr = strjoin(message, '\n');
    else
        messageStr = message;
    end
    fprintf('%s\n', messageStr);
end

%{
OLD CODE:

%}
