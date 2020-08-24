function [toQuit, toRedo, eventInfoOut, eventClassOut, isCheckedOut] = ...
                minEASE_examine_gapfree_events (eventInfoLast, eventClassLast, isCheckedLast, dataRaw, dataLowpass, tVec, outputDirectory, outputLabel, directionPsc, varargin)
%% Examine detected events and allow additions and deletions
% Usage: [toQuit, toRedo, eventInfoOut, eventClassOut, isCheckedOut] = ...
%               minEASE_examine_gapfree_events (eventInfoLast, eventClassLast, isCheckedLast, dataRaw, dataLowpass, tVec, outputDirectory, outputLabel, directionPsc, varargin)
% Explanation: 
%       TODO
% Outputs:
%       TODO
% Side Effects:
%       TODO
% Arguments:    
%       TODO
%       varargin    - 'NoiseLevel': Noise level for the signal to compare with
%                   must be a positive scalar
%                   default == root-mean-square level of Gaussian part of data
%                   - 'ZoomWindowMs': window for examining events (ms)
%                   must be a numeric positive scalar
%                   default == 50 ms
%                   - 'ToPrompt': whether to prompt for decisions
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   
%
% Requires:
%       cd/minEASE_gui_examine_events.m
%       cd/compute_rms_Gaussian.m
%
% Used by:
%       cd/minEASE.m

% File History:
% ---------- Created by Koji Takahashi & Mark P Beenhakker
% 2017-06-06 AL - Renamed goldfinger_minis.m -> minEASE_examine_gapfree_events.m
% 2017-06-12 AL - Passed directionPsc to GUI for adjust_peaks
% 2017-07-24 AL - Moved code to dlmwrite_with_header.m
% 2017-07-24 AL - Instead of using Auto versions, now uses Last versions
%                   of eventInfo, eventClass, isChecked
% 2017-07-25 AL - Now saves sweep labels too
% 2017-10-15 AL - Moved saving events to minEASE.m
% 2018-01-28 AL - Added isdeployed
% 2018-02-13 AL - Added NoPrompts
% 2018-02-16 AL - Now passes noiseLevel to minEASE_gui_examine_events.m
% 2018-02-23 AL - Now checks whether GUI is closed every 0.1 seconds
% 2018-02-27 AL - Changed NoPrompts to toPrompt
% 2018-08-03 AL - Renamed sweepLabel -> outputLabel
%

% Shared with minEASE_gui_examine_events.m
global eventInfo eventClass isChecked                   % used
global finalCommand                                     % used

%% Hard-coded parameters
PAUSE_TIME = 0.1;

%% Default parameters for TODO
zoomWindowMsDefault = 50;           % default window for examining events (ms)
toPromptDefault = true;             % whether to prompt for decisions by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 7
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NoiseLevel', [], ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ZoomWindowMs', zoomWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ToPrompt', toPromptDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
noiseLevel   = iP.Results.NoiseLevel;
zoomWindowMs = iP.Results.ZoomWindowMs;
toPrompt     = iP.Results.ToPrompt;

% If no noiseLevel provided, use the root mean square of the 
%   Gaussian part of the original data trace as the default
%   Note: this should not happen in minEASE
if isempty(noiseLevel)
    noiseLevel = compute_rms_Gaussian(data);            % default Gaussian noise level
end

%% Initialize outputs
toQuit = false;
toRedo = true;
eventInfo = eventInfoLast;
eventClass = eventClassLast;
isChecked = isCheckedLast;

%% Create and show GUI
%% TODO: Pass in other optional arguments
[hfig.gui] = minEASE_gui_examine_events(dataRaw, dataLowpass, tVec, ...
                                outputLabel, directionPsc, ...
                                'NoCloseReq', ~toPrompt, ...
                                'NoiseLevel', noiseLevel);

% Wait for user to close GUI before continuing
while ishandle(hfig.gui)    % while GUI is still there
    % Pause for PAUSE_TIME seconds
    pause(PAUSE_TIME);          
end

% Upon clicking "DONE!" (closes GUI), radiobutton tag is returned
switch finalCommand
case 'NextTrace'
    toRedo = false;
case 'StartOver'
    %toRedo = true;
case 'Quit'
    toQuit = true;
end

%% TODO: Don't close GUI but replot instead
% Make sure all other figures are closed
close all force

% Return new eventInfo, eventClass & isChecked
eventInfoOut = eventInfo;
eventClassOut = eventClass;
isCheckedOut = isChecked;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
