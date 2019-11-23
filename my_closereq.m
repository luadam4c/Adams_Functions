function [structure] = my_closereq (hObject, varargin)
%% Close request function that displays a question dialog box
% Usage: [structure] = my_closereq (hObject, varargin)
%   hObject         - handle of object
%   callbackData    - (opt) %%% TODO: examine
%   defaultOption   - (opt) the default option
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'Yes'       - close the object
%                       'No'        - don't close the object
%                   default == 'No'
%   objectName      - (opt) the name of the object to be closed
%                   must be a string scalar or a character vector
%
% Used by:
%       cd/m3ha_optimizergui_4compgabab.m
%       /home/Matlab/minEASE/gui_examine_events.m
%
% 2016-04-25 Adapted from Figure properties documentation
% 2016-06-16 Added input parser scheme and added defaultOption
% 

%% Default values for optional arguments
callbackDataDefault = [];
objectNameDefault = 'this figure';
defaultOptionDefault = 'No';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'my_closereq';

% Add required inputs to an Input Parser
addRequired(iP, 'hObject', @ishandle);      % handle of object

% Add optional inputs to the Input Parser
addOptional(iP, 'callbackData', callbackDataDefault);
addOptional(iP, 'defaultOption', defaultOptionDefault, ...
    @(x) any(validatestring(x, {'Yes', 'No'})));
addOptional(iP, 'objectName', objectNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Read from the Input Parser
parse(iP, hObject, varargin{:});
callbackData = iP.Results.callbackData;
defaultOption = validatestring(iP.Results.defaultOption, {'Yes', 'No'});
objectName = iP.Results.objectName;

%% Open question dialog box
selection = questdlg(sprintf('Are you SURE you want to close %s?', objectName), ...
            'Close Request Function', 'Yes', 'No', defaultOption);

%% Close hObject or not
switch selection
case 'Yes'
    delete(hObject);
case 'No'
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
